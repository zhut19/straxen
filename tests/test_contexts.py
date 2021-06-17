"""For all of the context, do a quick check to see that we are able to search
a field (i.e. can build the dependencies in the context correctly)
See issue #233 and PR #236"""
from straxen.contexts import *
import tempfile
import os


def import_wfsim():
    """Check if we can test the wfsim-related contexts. This is not the case
    for e.g. travis checks as we don't install wfsim by default."""
    try:
        import wfsim
        return True
    except ImportError:
        # We cannot test the wfsim as it is not installed.
        return False

##
# XENONnT
##


def test_xenonnt_online():
    st = xenonnt_online(_database_init=False, use_rucio=False)
    st.search_field('time')


def test_xennonnt():
    st = xenonnt(_database_init=False, use_rucio=False)
    st.search_field('time')


def test_xennonnt_latest(cmt_version='latest'):
    if straxen.utilix_is_configured():
        st = xenonnt(cmt_version, _database_init=False, use_rucio=False)
        st.search_field('time')


def test_xenonnt_led():
    st = xenonnt_led(_database_init=False, use_rucio=False)
    st.search_field('time')


def test_nt_is_nt_online():
    if not straxen.utilix_is_configured():
        # Cannot contact CMT without the database
        return
    # Test that nT and nT online are the same
    st_online = xenonnt_online(_database_init=False, use_rucio=False)

    st = xenonnt(_database_init=False, use_rucio=False)
    for plugin in st._plugin_class_registry.keys():
        print(f'Checking {plugin}')
        nt_key = st.key_for('0', plugin)
        nt_online_key = st_online.key_for('0', plugin)
        assert str(nt_key) == str(nt_online_key)


def test_xenonnt_simulation():
    if import_wfsim():
        st = xenonnt_simulation()
        st.search_field('time')

##
# XENON1T
##


def test_xenon1t_dali():
    st = xenon1t_dali()
    st.search_field('time')


def test_demo():
    """
    Test the demo context. Since we download the folder to the current
        working directory, make sure we are in a tempfolder where we
        can write the data to
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            print("Temporary directory is ", temp_dir)
            os.chdir(temp_dir)
            st = demo()
            st.search_field('time')
        # On windows, you cannot delete the current process'
        # working directory, so we have to chdir out first.
        finally:
            os.chdir('..')


def test_fake_daq():
    st = fake_daq()
    st.search_field('time')


def test_xenon1t_led():
    st = xenon1t_led()
    st.search_field('time')


def test_xenon1t_simulation():
    if import_wfsim():
        st = xenon1t_simulation()
        st.search_field('time')
