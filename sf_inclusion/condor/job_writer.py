# coding: utf-8

_all_ = [ "JobWriter" ]

import os
from functools import wraps

class JobWriter:
    """
    Help writing shell and condor job files.
    Writes one file at a time.
    """
    def __init__(self):
        self.filenames = []
        self.exts = ('sh', 'condor')
        self.endl = '\n'
        self.par = lambda s : s + self.endl

    def add_string(self, string):
        with open(self.filenames[-1], 'a') as self.f:
            self.f.write( self.par(string) )

    @staticmethod
    def define_output(localdir, data_folders, tag, names=None):
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
        if names is None:
            names = data_folders
        if not isinstance(names, (tuple,list)):
            names = [ names ]            
        if len(data_folders) > 1:
            assert len(data_folders) == len(names)

        # ensure the length of the two lists is the same for the `zip` that follows
        if len(data_folders) == 1 and len(names) > 1:
            data_folders = len(names)*data_folders

        job_d, check_d = ([] for _ in range(2))
        for dataf in data_folders:
            job_d.append( os.path.join(base_d, 'submission', dataf) )
            mkdir(job_d[-1])

            check_d.append( os.path.join(base_d, 'outputs', dataf) )
            mkdir(check_d[-1])

        job_f, subm_f, check_f = ([] for _ in range(3))
        for jd, cd, name in zip(job_d,check_d,names):
            mkdir(cd)
            job_f.append( os.path.join(jd, name + '.sh') )
            subm_f.append( os.path.join(jd, name + '.condor') )
            
            check_name = '{name}_C$(Cluster)P$(Process).o'
            check_name = check_name.format(name=name)
            check_f.append( os.path.join(cd, check_name) )
                
        return job_f, subm_f, check_f

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

    def write_condor(self, filename, executable, outfile, queue, machine='llrt3condor'):
        self.filenames.append( filename )
        m = ( self.par('Universe = vanilla') +
              self.par('Executable = {}'.format(executable)) +
              self.par('input = /dev/null') +
              self.par('output = {}'.format(outfile)) +
              self.par('error  = {}'.format(outfile.replace('.o', '.e'))) +
              self.par('getenv = true') +
              self.condor_specific_content(queue=queue, machine=machine) +
              self.endl )
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
        m = ( self.par('T3Queue = {}'.format(queue)) +
              self.par('WNTag=el7') +
              self.par('+SingularityCmd = ""') )
        if machine == 'llrt3condor':
            m += self.endl + 'include : /opt/exp_soft/cms/t3/t3queue |'
            # m += ( self.par('UNIX_GROUP = cms') +
            #        self.par('if !defined T3Queue') +
            #        self.par('  T3Queue=short') +
            #        self.par('endif') +

            #        self.par('accounting_group = $(UNIX_GROUP)') +

            #        self.par('concurrency_limits_expr = strcat(T3Queue,":",RequestCpus," ",AcctGroupUser,":",RequestCpus)') +

            #        self.par('+T3Queue="$(T3Queue)"') +
            #        self.par('+T3Group="$(UNIX_GROUP)"') +
            #        self.par('+WNTag="$(WNTag)"') +
            #        self.par('T3Submit=true') )
                  
        elif machine == 'llrt3condor7':
            m += self.endl + 'include : /opt/exp_soft/cms/t3_tst/t3queue |'
        else:
            raise ValueError('Machine {} is not supported.'.format(machine))
        m += self.endl
        return m
