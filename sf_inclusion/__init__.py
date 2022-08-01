# from luigi_conf import luigi_cfg
import condor
import scripts
import luigi_conf

os.environ['LUIGI_CONFIG_PATH'] = os.path.join(os.environ['PWD'], 'sf_inclusion/luigi_conf/luigi.cfg')
assert os.path.exists(os.environ['LUIGI_CONFIG_PATH'])

