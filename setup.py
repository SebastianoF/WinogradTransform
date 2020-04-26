from pip._internal.network.session import PipSession
from pip._internal.req import parse_requirements
from setuptools import find_packages, setup

base_requirements = [
    str(requirement).split()[0]
    for requirement in parse_requirements("requirements.txt", session=PipSession())
]


setup(
    name="Winograd",
    version="v0.0.0",
    description="Winograd transform Experiments",
    author="sebastiano ferraris",
    author_email="sebastiano.ferraris@gmail.com",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    license="None",
    python_requires=">=3.7",
    url="https://github.com/SebastianoF/WinogradTransform",
    packages=find_packages(include=["winograd"], exclude=["docs", "examples", "tests"]),
    install_requires=base_requirements,
)
