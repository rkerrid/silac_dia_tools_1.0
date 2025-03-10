from setuptools import setup, find_packages

setup(
    name='silac_dia_tools',
    version='0.2.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas>=1.0.0',  # specify the minimum version you require
        'matplotlib>=3.0.0',
        'numpy>=1.18.0',
        'jupyterlab>=3.0.0',
        'spyder',
        'icecream',
        'seaborn',
        'fpdf',
        'pytest',
        'tqdm',
        'pandastable',       
        'alphapept @ git+https://github.com/MannLabs/alphapept.git',  # specify the branch, tag, or commit
        "directlfq @ git+https://github.com/MannLabs/directlfq.git",# specify the branch, tag, or commit
        
    ],
)


# need to install pandastable conda install -c conda-forge pandastable
