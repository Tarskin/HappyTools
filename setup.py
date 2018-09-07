from setuptools import setup, find_packages
setup(
    name="HappyTools",
    version="0.0.2",
    packages=find_packages(),
    install_requires=[
        'scipy>=1.1.0',
        'matplotlib>=2.2.3',
        'numpy>=1.15.1'
    ],

    # metadata to display on PyPI
    author="Bas Cornelis Jansen",
    author_email="bas.c.jansen@gmail.com",
    description="HappyTools (Python3) software package for chromatography",
    license="Apache2",
    keywords="HappyTools",
    url="https://github.com/Tarskin/HappyTools",
    project_urls={
        "Bug Tracker": "https://github.com/Tarskin/HappyTools/issues",
        "Source Code": "https://github.com/Tarskin/HappyTools",
    }
)
