import setuptools

PKG_NAME = "moldo"
PKG_DIR = PKG_NAME.replace('-', '_')

with open(f"./src/{PKG_NAME}/version.txt") as f:
    version = f.read()
    
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

pks = [PKG_DIR]
setuptools.setup(
    name=PKG_NAME,
    version=version,
    author="Tony, Avery, Patrik",
    author_email="contacttonyliu@gmail.com",
    description="mini molecular docking workflow",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=f"https://github.com/Tony-xy-Liu/{PKG_NAME}",
    project_urls={
        "Bug Tracker": f"https://github.com/Tony-xy-Liu/{PKG_NAME}/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    # packages=setuptools.find_packages(where="src"),
    packages=pks,
    # package_data={
    #     # "":["*.txt"],
    #     # "package-name": ["*.txt"],
    #     # "test_package": ["res/*.txt"],
    # },
    entry_points={
        'console_scripts': [
            f'moldo = {PKG_DIR}:main',
        ]
    },
    python_requires=">=3.9",
)
