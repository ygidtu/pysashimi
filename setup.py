from setuptools import setup, find_packages

setup(
    name='pySashimi',
    author='ChenlinLab',
    author_email='ygidtu@gmail.com',
    version='1.5.0',
    long_description=__doc__,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    maintainer='ygidtu',
    maintainer_email='ygidtu@gmail.com',
    url='',
    entry_points={
        'console_scripts':
            [
                'pysashimi = cli.cli:cli'
            ]
    },
    data_files=[(".", ['settings.ini'])],
    # scripts=['main.py'],
    install_requires=[
        "click",
        "cycler",
        "et-xmlfile",
        "filetype",
        "jdcal",
        "kiwisolver",
        "matplotlib",
        "numpy",
        "openpyxl",
        "pyparsing",
        "pysam",
        "python-dateutil",
        "six",
        "tqdm",
    ],
)