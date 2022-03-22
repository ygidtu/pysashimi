import os

from setuptools import setup, find_packages


def locate_packages():
    __dir__ = os.path.dirname(os.path.abspath(__file__))
    pipfile = os.path.join(__dir__, "Pipfile")
    packages = []
    if os.path.exists(pipfile):
        kept = False
        with open(pipfile) as r:
            for line in r:
                line = line.strip()

                if line == "[packages]":
                    kept = True
                    continue
                elif line.startswith("[") and line.endswith("]"):
                    kept = False

                if kept:
                    if line:
                        name, version = line.split(" = ")
                        version = version.strip('"')
                        if version == "*":
                            version = ""
                        packages.append("%s %s" % (name, version))
    return packages


setup(
    name='pySashimi',
    author='ygidtu',
    author_email='ygidtu@gmail.com',
    version='1.6.0',
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
                'pysashimi = cli.cli:plot'
            ]
    },
    python_requires='>=3.6',
    data_files=[(".", ['settings.ini'])],
    install_requires=locate_packages(),
)
