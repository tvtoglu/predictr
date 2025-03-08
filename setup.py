from setuptools import setup

with open("README.md", "r") as fh:
          long_description = fh.read()

setup(
      name='predictr',
      version='0.1.28',
      description='Life Data Analysis for Reliability Engineers - Weibull Analysis, Detailed Plots, Compute Statistics',
      author='Tamer Tevetoglu',
      author_email="predictr@outlook.com",
      url="https://tvtoglu.github.io/predictr/",
      project_urls={
        'Documentation': "https://tvtoglu.github.io/predictr/classes/",
        'GitHub': "https://github.com/tvtoglu/predictr",
        'Github Discussions': "https://github.com/tvtoglu/predictr/discussions"},
      py_modules=["predictr"],
      package_dir={'': 'src'},
      long_description=long_description,
      long_description_content_type='text/markdown',
      keywords = 'reliability, weibull, bias, life data analysis, engineering, confidence, bootstrap, monte-carlo, fisher bounds, likelihood ratio, unreliability, survival analysis, lifelines, testing',
      classifiers=[
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
          "Programming Language :: Python :: 3.10",
          "Programming Language :: Python :: 3.11",
          "Programming Language :: Python :: 3.12",
          "Programming Language :: Python :: 3.13",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",],
      install_requires= [
          "numpy >= 1.16.0",
          "scipy >= 1.3.0",
          "pandas >= 1.0.0",
          "matplotlib >= 2.2.0"
          ]
      )
