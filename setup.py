import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vk2gpz.geom",
    version="1.0.3",
    author="vk2gpz",
    author_email="vk2gpz@gmail.com",
    license='vk2gpz',
    description="A collection of geometry related modules.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vk2gpz/vk_pygeom",
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_namespace_packages(include=['vk2gpz.*']),
    python_requires=">=3.9",
)
