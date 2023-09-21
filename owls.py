import os
import solara

@solara.component
def Page():
    import numpy as np

    import astropy.units as u
    from astropy.utils.data import download_file
    from astropy.nddata import StdDevUncertainty

    from specutils import SpectrumList, SpectrumCollection

    from glue.core.roi import XRangeROI
    from glue.core.edit_subset_mode import NewMode

    import ipyvuetify as v
    import solara

    from jdaviz.core.helpers import ConfigHelper
    import jdaviz
    from jdaviz.app import custom_components
    from specutils import Spectrum1D

    import json
    import ipyvue
    import ipysplitpanes
    import ipygoldenlayout

    ipysplitpanes.SplitPanes()
    ipygoldenlayout.GoldenLayout()
    for name, path in custom_components.items():
        ipyvue.register_component_from_file(None, name,
                                            os.path.join(os.path.dirname(jdaviz.__file__), path))

    ipyvue.register_component_from_file('g-viewer-tab', "container.vue", jdaviz.__file__)

    all_owls_links_path = 'all_owls_links.json'

    links = json.load(open(all_owls_links_path, 'r'))

    class OWLSviz(ConfigHelper):
        """OWLS-Specviz Helper class."""

        _default_configuration = "./owls.yaml"

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._default_spectrum_viewer_reference_name = None

        def load_data(self, data, data_label=None, format=None, show_in_viewer=True,
                      concat_by_file=False):
            """
            Load data into OWLSviz.

            Parameters
            ----------
            data : str, `~specutils.Spectrum1D`, or `~specutils.SpectrumList`
                Spectrum1D, SpectrumList, or path to compatible data file.
            data_label : str
                The Glue data label found in the ``DataCollection``.
            format : str
                Loader format specification used to indicate data format in
                `~specutils.Spectrum1D.read` io method.
            show_in_viewer : bool
                Show data in viewer(s).
            concat_by_file : bool
                If True and there is more than one available extension, concatenate
                the extensions within each spectrum file passed to the parser and
                add a concatenated spectrum to the data collection.
            """
            super().load_data(data,
                              parser_reference='specviz-spectrum1d-parser',
                              data_label=data_label,
                              format=format,
                              show_in_viewer=show_in_viewer,
                              concat_by_file=concat_by_file)

        def get_data(self, data_label=None, spectral_subset=None, cls=None,
                     use_display_units=False, **kwargs):
            """
            Returns data with name equal to data_label of type cls with subsets applied from
            spectral_subset.

            Parameters
            ----------
            data_label : str, optional
                Provide a label to retrieve a specific data set from data_collection.
            spectral_subset : str, optional
                Spectral subset applied to data.
            cls : `~specutils.Spectrum1D`, optional
                The type that data will be returned as.
            use_display_units: bool, optional
                Whether to convert to the display units defined in the <unit-conversion> plugin.

            Returns
            -------
            data : cls
                Data is returned as type cls with subsets applied.

            """
            spatial_subset = kwargs.pop("spatial_subset", None)
            function = kwargs.pop("function", None)
            if len(kwargs) > 0:
                raise ValueError(f'kwargs {[x for x in kwargs.keys()]} are not valid')

            if cls is None:
                cls = Spectrum1D
            elif spatial_subset or function:
                raise ValueError('kwargs spatial subset and function are not valid in specviz')
            else:
                spatial_subset = None
                function = None

            return self._get_data(data_label=data_label, spatial_subset=spatial_subset,
                                  spectral_subset=spectral_subset, function=function,
                                  cls=cls, use_display_units=use_display_units)

    s = OWLSviz()
    # s.plugins['Model Fitting']._obj.dataset._viewers = ['h', 'k']

    h_centroid = 3968.4673 * u.Angstrom
    k_centroid = 3933.6614 * u.Angstrom
    roi_half_width = 3 * u.Angstrom

    rois = []
    for line in [k_centroid, h_centroid]:
        rois.append(
            XRangeROI(
                (line - roi_half_width).to(u.AA).value,
                (line + roi_half_width).to(u.AA).value
            )
        )

    def add_spectrum_from_url(viz_helper, data_url):
        all_orders = SpectrumCollection.read(
            download_file(data_url, cache=True)
        )
        target_name = all_orders.meta['header']['OBJNAME']        

        # Delete existing subsets
        for subset_grp in viz_helper.app.data_collection.subset_groups:
            viz_helper.app.data_collection.remove_subset_group(subset_grp)

        for ref, order, roi in zip(['h', 'k'], [90, 89], rois):
            # viz_helper.app.state.viewer_icons[ref] = f"{order}"
            spectrum = all_orders[order]
            spectrum.uncertainty = StdDevUncertainty(np.sqrt(spectrum.flux.value))
            data_label = f"{target_name}, Order {order}"
            s._default_spectrum_viewer_reference_name = ref
            viz_helper.load_data(spectrum, data_label=data_label, show_in_viewer=ref)

            viz_helper.app.add_data_to_viewer(ref, data_label, clear_other_data=True)    
            spectrum_viewer = viz_helper.app.get_viewer(ref)

            # Set the active edit_subset_mode to NewMode to be able to add multiple subregions
            spectrum_viewer.session.edit_subset_mode._mode = NewMode
            spectrum_viewer.apply_roi(roi)

    def on_click(widget, event, data):
        add_spectrum_from_url(s, links[data]['url'])
    
    
    with solara.Column():
        a = v.Select(
            v_model='HAT-P-11 @ 2021-09-11 07:27',
            label='Observation',
            items=sorted(links.keys()))

        a.on_event('change', on_click)

        add_spectrum_from_url(s, links[a.v_model]['url'])

        title = v.Html(
            tag='h1',
            children=['Olin Wilson Legacy Survey (OWLS) Spectrum Visualizer']
        )


        subtitle = v.Text(
            children=['Powered by jdaviz']
        )


        content_main = v.Layout(
            column=True,
            _metadata={'mount_id': 'content-main'},
            children=[title, subtitle, a, s.app]
        )

        display(content_main)
