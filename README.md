# jdaviz-owls
OWLS spectra interactively visualized via render.com

### Setting up for render.com

* Create an account (recommend signing in with GitHub)
* Click the colorful `New +` button and select `Web Service`
* Select your GitHub repository containing your solara app
* Under the `Environment` tab, set an environment variable with the key `PYTHON_VERSION` to `3.9.10`.
* Under the `Settings` tab, set `Start Command` to `solara run owls.py --port=$PORT --no-open --host=0.0.0.0`
