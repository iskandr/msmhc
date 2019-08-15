# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import (absolute_import,)

import os
import logging
import re

from setuptools import setup, find_packages

readme_dir = os.path.dirname(__file__)
readme_path = os.path.join(readme_dir, 'README.md')

try:
    with open(readme_path, 'r') as f:
        readme_markdown = f.read()
except:
    logging.warning("Failed to load %s" % readme_path)
    readme_markdown = ""


with open('msmhc/__init__.py', 'r') as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        f.read(),
        re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")

if __name__ == '__main__':
    setup(
        name='msmhc',
        version=version,
        description="Looking for dark matter of the MHC ligandome in mass spectra",
        author="Alex Rubinsteyn",
        author_email="alex.rubinsteyn@gmail.com",
        url="https://github.com/iskandr/msmhc",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            'pyensembl>=1.5.0',
            'varcode',
            'progressbar2',
            'pandas'
        ],
        long_description=readme_markdown,
        long_description_content_type='text/markdown',
        packages=find_packages(),
        entry_points={
            'console_scripts': [
                'msmhc-generate=msmhc.generate_cli:run',
                'msmhc-fdr=msmhc.fdr_cli:run'
            ]
        }
    )
