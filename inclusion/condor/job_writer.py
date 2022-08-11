# coding: utf-8

_all_ = [ 'JobWriter' ]

import os
from functools import wraps

import os
import sys
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, parent_dir)

import inclusion
from inclusion import config

class JobWriter:
    """
    Help writing shell and condor job files.
    Writes one file at a time.
    """
    def __init__(self):
        self.filenames = []
        self.exts = ('sh', 'condor')
        self.endl = '\n'

    def add_string(self, string):
        with open(self.filenames[-1], 'a') as self.f:
            self.f.write( string + self.endl )

    @staticmethod
    def define_output(localdir, data_folders, tag, names=''):
        """
        Defines where the shell and condor job files, and the HTCondor outputs
        will be stored.
        - When `data_folders` has more than one element, `names` will be matched in
          terms of length, so that each name will correspond to a different folder. This
          has a consequence on the HTCondor output files only.
        """
        base_d = os.path.join(localdir, 'jobs', tag)
        mkdir = lambda d : os.system('mkdir -p {}'.format(d))

        # type checks
        if not isinstance(data_folders, (tuple,list)):
            data_folders = [ data_folders ]
        if not isinstance(names, (tuple,list)):
            names = [ names ]            
        if len(data_folders) != len(names):
            if len(names) == 1:
                names = len(data_folders)*names
            else:
                raise ValueError('You got the total number of foldernames wrong.')
        

        # ensure the length of the two lists is the same for the `zip` that follows
        if len(data_folders) == 1 and len(names) > 1:
            data_folders = len(names)*data_folders

        job_d, out_d = ([] for _ in range(2))
        for dataf in data_folders:
            job_d.append(os.path.join(base_d,
                                      config.folders['subm'], dataf))
            mkdir(job_d[-1])
            out_d.append( os.path.join(base_d,
                                       config.folders['outs'], dataf) )
            mkdir(out_d[-1])

        job_f, subm_f, out_f, log_f = ([] for _ in range(4))
        for jd, cd, name in zip(job_d,out_d,names):
            mkdir(cd)
            job_f.append( os.path.join(jd, 'job{}.sh'.format(name)) )
            subm_f.append( os.path.join(jd, 'job{}.condor'.format(name)) )

            base_name = 'C$(Cluster)_P$(Process)'
            out_name = '{}.out'.format(base_name)
            log_name = '{}.log'.format(base_name)
            out_f.append( os.path.join(cd, out_name) )
            log_f.append( os.path.join(cd, log_name) )
                
        return job_f, subm_f, out_f, log_f

    @staticmethod
    def define_dag_output(localdir, tag, name):
        assert '.' not in name
        mkdir = lambda d : os.system('mkdir -p {}'.format(d))
        subm_d = os.path.join(localdir, 'jobs', tag, 'outputs')
        mkdir(subm_d)
        subm_d = os.path.join(subm_d, 'CondorDAG')
        mkdir(subm_d)
        return os.path.join(subm_d, name + '.dag')

    def extension_exception(self):
        mes = 'Wrong extension. Supported ones are:'
        for ext in self.exts:
            mes += '  - {}'.format(ext)
        raise ValueError(mes)
        
    def check_mode_(self, mode):
        if mode not in self.exts:
            raise ValueError('The mode {} is not supported.'.format(mode))

    def write_condor(self, filename, shell_exec, real_exec, outfile, logfile,
                     queue, machine='llrt3condor'):
        self.filenames.append(filename)
        batch_name = os.path.splitext(os.path.basename(executable))[0]
        m = self.endl.join(('Universe = vanilla',
                            'Executable = {}'.format(shell_exec),
                            'input = {}'.format(real_exec),
                            'output = {}'.format(outfile),
                            'error = {}'.format(outfile.replace('.out', '.err')),
                            'log = {}'.format(logfile),
                            'getenv = true',ful
                            '+JobBatchName="{}"'.format(batch_name),
                            'should_transfer_files = YES',
                            self.condor_specific_content(queue=queue, machine=machine)))
        m += self.endl
        with open(filename, 'w') as self.f:
            self.f.write(m)
        os.system('chmod u+rwx '+ filename)

    def write_queue(self, qvars=(), qlines=[]):
        """
        Works for any number variables in the queue.
        It is up to the user to guarantee compatibility between queue variables and lines.
        """
        extension = self.filenames[-1].split('.')[-1]
        if extension != self.exts[1]:
            self.extension_exception()
            
        with open(self.filenames[-1], 'a') as self.f:
            if len(qvars) > 0:
                argstr = 'Arguments = '
                for qvar in qvars:
                    argstr += '$(' + qvar + ') '
                argstr += self.endl
                self.f.write(argstr)
                qstr = 'queue '
                for qvar in qvars:
                    qstr += ('' if qvar == qvars[0] else ',') + qvar
                qstr += ' from ('
                self.f.write(qstr + self.endl)
                for ql in qlines:
                    self.f.write(ql + self.endl)
                self.f.write(')' + self.endl)
            else:
                self.f.write('queue' + self.endl)

    def write_shell(self, filename, command, localdir):
        self.filenames.append( filename )
        m = ( '#!/bin/bash' +
              self.endl + 'export X509_USER_PROXY=~/.t3/proxy.cert' +
              self.endl + 'export EXTRA_CLING_ARGS=-O2' +
              self.endl + 'source /cvmfs/cms.cern.ch/cmsset_default.sh' +
              self.endl + 'cd {}/'.format(localdir) +
              self.endl + 'eval `scramv1 runtime -sh`' +
              self.endl + command +
              self.endl )
        with open(filename, 'w') as self.f:
            self.f.write(m)
        os.system('chmod u+rwx '+ filename)

    def condor_specific_content(self, queue, machine='llrt3condor'):
        m = self.endl.join(('T3Queue = {}'.format(queue),
                            'WNTag=el7',
                            '+SingularityCmd = ""'))
        if machine == 'llrt3condor':
            m += self.endl + 'include : /opt/exp_soft/cms/t3/t3queue |'
                  
        elif machine == 'llrt3condor7':
            m += self.endl + 'include : /opt/exp_soft/cms/t3_tst/t3queue |'
        else:
            raise ValueError('Machine {} is not supported.'.format(machine))
        m += self.endl
        return m
