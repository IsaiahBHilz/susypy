import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name='susypy',
    version='1.0.0',
    author='Saul & Isaiah B. Hilsenrath',
    author_email='ihilsenr@umd.edu',
    description='A symbolic algebra system for supersymmetry calculations.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/IsaiahBHilz/susypy',
    packages=setuptools.find_packages(),
    install_requires=['sympy', 'cadabra2'],
    license='GPLv3',
    license_files = ('COPYING.txt',),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3 :: Only',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.8'
)