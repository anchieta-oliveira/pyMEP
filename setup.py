from setuptools import setup, find_packages

def readme():
    with open("README.md", encoding="utf-8") as f:
        return f.read()

setup(
    name="pymep",
    version="0.0.1",
    description="LabMonitor, a Python application designed to simplify the management of computing resources in decentralized networks of Linux machines.",
    long_description=readme(),
    long_description_content_type="text/markdown",
    author_email="anchieta.oliveira@biof.ufrj.br",
    url="https://github.com/anchieta-oliveira/pyMEP",
    license="MIT",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19.2",
        "streamlit>=1.41.1",
        "sphinx>=8.1.3",
        "plotly>=5.24.1",
    ],
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
