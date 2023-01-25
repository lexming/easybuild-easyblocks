##
# Copyright 2021 Vrije Universiteit Brussel
#
# This file is part of EasyBuild,
# originally created by the HPC team of Ghent University (http://ugent.be/hpc/en),
# with support of Ghent University (http://ugent.be/hpc),
# the Flemish Supercomputer Centre (VSC) (https://www.vscentrum.be),
# Flemish Research Foundation (FWO) (http://www.fwo.be/en)
# and the Department of Economy, Science and Innovation (EWI) (http://www.ewi-vlaanderen.be/en).
#
# https://github.com/easybuilders/easybuild
#
# EasyBuild is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation v2.
#
# EasyBuild is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with EasyBuild.  If not, see <http://www.gnu.org/licenses/>.
##
"""
EasyBuild support for building and installing SOWFA, implemented as an easyblock

@author: Alex Domingo (Vrije Universiteit Brussel)
"""
import os

from easybuild.easyblocks.generic.cmdcp import MakeCp
from easybuild.framework.easyconfig import CUSTOM
from easybuild.tools.filetools import find_glob_pattern
from easybuild.tools.run import run_cmd

class EB_SOWFA(MakeCp):
    """Support for building and installing SOWFA."""

    @staticmethod
    def extra_options():
        """Change default values of options"""
        extra = MakeCp.extra_options()
        # files_to_copy is not mandatory here
        extra['files_to_copy'][2] = CUSTOM
        return extra

    def build_step(self):
        """Build SOWFA with wmake."""
        sowfa_dir_pattern = os.path.join(self.builddir, '%s-*' % self.name)
        sowfa_dir = find_glob_pattern(sowfa_dir_pattern)

        build_env = {
            'SOWFA_DIR': sowfa_dir,
            'WM_PROJECT_USER_DIR': sowfa_dir,
            'WM_NCOMPPROCS': self.cfg['parallel'],
        }

        cmd = 'source $FOAM_BASH && '
        cmd += ' && '.join(['export %s="%s"' % (k, build_env[k]) for k in build_env])
        cmd += ' && ./Allwclean && ./Allwmake'

        return run_cmd(cmd, log_all=True, simple=True, log_output=True)

    def install_step(self):
        """Copy files to installation directory."""
        self.cfg['files_to_copy'] = ['applications', 'exampleCases', 'lib', 'tools', 'README.SOWFA']

        super(EB_SOWFA, self).install_step()

        # catch $WM_OPTIONS from library path
        sowfa_libdir_pattern = os.path.join(self.installdir, 'lib', '*')
        self.wm_options = os.path.basename(find_glob_pattern(sowfa_libdir_pattern))

    def sanity_check_step(self):
        """Custom sanity check for SOWFA."""
        sowfa_bin = ['ABLSolver', 'ABLTerrainSolver', 'pisoFoamTurbine.ALM', 'setFieldsABL',
                     'turbineTestHarness.ALM', 'windPlantSolver.ALM']
        app_bin_path = os.path.join('applications', 'bin', self.wm_options)
        sowfa_bin_paths = [os.path.join(app_bin_path, x) for x in sowfa_bin]

        custom_paths = {
            'files': ['README.SOWFA'] + sowfa_bin_paths,
            'dirs': [os.path.join('lib', self.wm_options), 'exampleCases', 'tools'],
        }

        super(EB_SOWFA, self).sanity_check_step(custom_paths=custom_paths)

    def make_module_req_guess(self):
        """Set the environment."""
        guesses = super(EB_SOWFA, self).make_module_req_guess()

        guesses['PATH'] = [os.path.join('applications', 'bin', self.wm_options)]
        guesses['LIBRARY_PATH'] = [os.path.join('lib', self.wm_options)]
        guesses['LD_LIBRARY_PATH'] = [os.path.join('lib', self.wm_options)]

        return guesses
