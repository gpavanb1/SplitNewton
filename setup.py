from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='SplitNewton',
      version='0.1.3',
      description='Split Newton Solver',
      url='https://github.com/gpavanb1/SplitNewton',
      author='gpavanb1',
      author_email='gpavanb@gmail.com',
      license='MIT',
      packages=['splitnewton'],
      install_requires=["numpy",
      "scipy"],
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=[
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
      ],
      keywords='newton python continuation armijo optimization pseudotransient splitting',
      project_urls={  # Optional
        'Bug Reports': 'https://github.com/gpavanb1/SplitNewton/issues',
        'Source': 'https://github.com/gpavanb1/SplitNewton/',
      },
      zip_safe=False)