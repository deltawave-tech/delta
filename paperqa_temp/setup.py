from setuptools import setup, find_packages

setup(
    name="paperqa",
    version="5.8.0",
    packages=find_packages(),
    install_requires=[
        "aiohttp",
        "anyio",
        "coredis",
        "fhaviary",
        "html2text",
        "httpx",
        "limits",
        "litellm",
        "numpy",
        "pybtex",
        "pydantic",
        "pydantic-settings",
        "PyMuPDF",
        "rich",
        "setuptools",
        "tantivy",
        "tenacity",
        "tiktoken"
    ],
)