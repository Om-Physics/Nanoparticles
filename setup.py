#!/usr/bin/env python3
"""
setup.py
========
Package configuration for the gold nanoparticle drug delivery
computational framework.

Install in development mode:
    pip install -e .

Install normally:
    pip install .
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [
        line.strip()
        for line in fh
        if line.strip() and not line.startswith("#")
    ]

setup(
    name="aunp-drug-delivery",
    version="1.0.0",
    author="Om Jha",
    author_email="om.physics7@gmail.com",
    description=(
        "Computational framework for gold nanoparticle drug delivery "
        "analysis in lung cancer — pharmacokinetics, biodistribution, "
        "toxicity, and interactive simulation."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Om-Physics/Nanoparticles",
    project_urls={
        "Bug Tracker": "https://github.com/Om-Physics/Nanoparticles/issues",
        "Documentation": "https://github.com/Om-Physics/Nanoparticles/tree/main/docs",
    },
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.8",
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Intended Audience :: Science/Research",
        "Development Status :: 5 - Production/Stable",
    ],
    entry_points={
        "console_scripts": [
            "generate-figures=generate_figures:main",
        ],
    },
    keywords=[
        "gold nanoparticles", "drug delivery", "lung cancer",
        "pharmacokinetics", "EPR effect", "computational modelling",
        "nanomedicine", "EGFR targeting",
    ],
)
