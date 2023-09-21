import os
import solara

@solara.component
def Page():
    import numpy as np

    import astropy.units as u
    from astropy.utils.data import download_file

    from specutils import SpectrumCollection, SpectralRegion

    import ipyvuetify as v

    import jdaviz
    from jdaviz import Specviz
    from jdaviz.app import custom_components
    from specutils import Spectrum1D

    import json
    import ipyvue
    import ipysplitpanes
    import ipygoldenlayout

    from specutils.fitting import fit_generic_continuum
    from specutils.manipulation import trapezoid_smooth
    from astropy.modeling.models import Chebyshev1D
    from astropy.stats import mad_std
    from scipy.stats import binned_statistic
    from specutils.manipulation import extract_region
    
    urls_loaded = set()
    def combine_orders(all_orders):
        normalized_orders = []
        continuum_orders = []
        for original_spectrum, exclude_regions in zip(
            [all_orders[89], all_orders[90]], 
            [SpectralRegion(3950*u.AA, 3975*u.AA), SpectralRegion(3915*u.AA, 3950*u.AA)]
        ):
            spectrum = trapezoid_smooth(original_spectrum, 201)
            x = spectrum.spectral_axis
            model = Chebyshev1D(6)
            g1_fit = fit_generic_continuum(spectrum, model=model, median_window=1, exclude_regions=exclude_regions)
            y_continuum_fitted = g1_fit(x)

            sort_lam = np.argsort(original_spectrum.wavelength)
            normalized_spectrum = Spectrum1D(flux=(original_spectrum.flux / y_continuum_fitted)[sort_lam], spectral_axis=original_spectrum.wavelength[sort_lam])
            normalized_orders.append(normalized_spectrum)
            continuum_orders.append(y_continuum_fitted[sort_lam])

        extracted_spectra = []

        for i in range(2):
            bs_std = binned_statistic(
                normalized_orders[i].wavelength.value,
                normalized_orders[i].flux.value,
                bins=100,
                statistic=lambda x: mad_std(x, ignore_nan=True)
            )

            bs_model = binned_statistic(
                normalized_orders[i].wavelength.value,
                continuum_orders[i],
                bins=100,
            )
            wl_bins = 0.5 * (bs_std.bin_edges[1:] + bs_std.bin_edges[:-1])

            p = np.polyfit(wl_bins - wl_bins.mean(), bs_std.statistic, 2, w=bs_model.statistic)
            fit = np.polyval(p, wl_bins - wl_bins.mean())

            min_std = fit.min()

            diff_sign = np.sign(fit - 2.5 * min_std)
            crossing_points = np.squeeze(np.argwhere(diff_sign[1:] != diff_sign[:-1]))

            high_snr_region = SpectralRegion(
                *(wl_bins[crossing_points] * normalized_orders[0].wavelength.unit)
            )
            extracted = extract_region(normalized_orders[i], high_snr_region)
            extracted_spectra.append(extracted)

        result = []
        flux_medians = []
        for i, order in enumerate(extracted_spectra):

            if i == 0:
                condition = order.wavelength.value < extracted_spectra[1-i].wavelength.value.max()
            else:
                condition = order.wavelength.value > extracted_spectra[1-i].wavelength.value.min()

            flux_median = np.nanmedian(order.flux.value[condition])
            flux_medians.append(flux_median)

            wl_span = order.wavelength.value[condition]

            normalization = 1
            if i == 1:
                normalization = flux_medians[0] / flux_median

            if i == 0:
                use_region = np.ones_like(order.flux.value).astype(bool)
            else: 
                use_region = np.logical_not(condition)

            #only use the non-overlapping wavelength region:
            result.append(
                np.column_stack([order.flux.value * normalization, order.wavelength.value])[use_region]
            )
        result = np.vstack(result)
        sort = np.argsort(result[:, 1])

        combined_orders = Spectrum1D(
            flux=result[sort, 0]*order.flux.unit, 
            spectral_axis=result[sort, 1]*order.wavelength.unit
        )
        return combined_orders
            
    
    ipysplitpanes.SplitPanes()
    ipygoldenlayout.GoldenLayout()
    for name, path in custom_components.items():
        ipyvue.register_component_from_file(None, name,
                                            os.path.join(os.path.dirname(jdaviz.__file__), path))

    ipyvue.register_component_from_file('g-viewer-tab', "container.vue", jdaviz.__file__)

    all_owls_links_path = 'all_owls_links.json'

    links = json.load(open(all_owls_links_path, 'r'))
            
    s = Specviz()
    s.app.template.template = s.app.template.template.replace(
        'calc(100% - 48px);', '800px'#'80vh !important;'
    )
    
    height = '800px'
    s.app.layout.height = height
    s.app.state.settings['context']['notebook']['max_height'] = height

    h_centroid = 3968.4673 * u.Angstrom
    k_centroid = 3933.6614 * u.Angstrom
    roi_half_width = 3 * u.Angstrom

    def add_spectrum_from_url(viz_helper, data_url):
        if data_url not in urls_loaded:
            all_orders = SpectrumCollection.read(
                download_file(data_url, cache=True)
            )
            target_name = all_orders.meta['header']['OBJNAME']        
            date = ':'.join(all_orders.meta['header']['DATE-OBS'].split(':')[:-1])

            data_label = f'{target_name} @ {date}'
            viz_helper.load_data(combine_orders(all_orders), data_label)
            spectrum_viewer = viz_helper.app.get_viewer(viz_helper._default_spectrum_viewer_reference_name)            
            spectrum_viewer.state.reset_limits()
            urls_loaded.add(data_url)
        
    def on_click(data):
        for link_key in data:
            add_spectrum_from_url(s, links[link_key]['url'])
    
    with solara.Column():
        solara.Markdown("# [Olin Wilson Legacy Survey (OWLS)](https://owls.readthedocs.io/)\n### Interactive spectrum visualizer")

        link_selected = solara.reactive([])
        solara.SelectMultiple("Observation", link_selected, list(links.keys()), on_value=on_click)
        
        display(s.app)

        solara.Markdown(
            f"CaII K and H = {k_centroid.value:.2f} and {h_centroid:.2f} (rest)\n\n"
            "Powered by [jdaviz](https://jdaviz.readthedocs.io/) and [solara](https://solara.dev/)"
        )
