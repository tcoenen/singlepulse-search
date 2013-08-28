from distutils.core import setup

setup(
    name='SSPS',
    version='0.2',
    packages=['ssps'],
    long_description=open('README.md').read(),
    scripts=['scripts/ssps_condense.py', 'scripts/ssps_grab.py'],
)
