##
# Copyright 2009-2013 Ghent University
#
# This file is part of EasyBuild,
# originally created by the HPC team of Ghent University (http://ugent.be/hpc/en),
# with support of Ghent University (http://ugent.be/hpc),
# the Flemish Supercomputer Centre (VSC) (https://vscentrum.be/nl/en),
# the Hercules foundation (http://www.herculesstichting.be/in_English)
# and the Department of Economy, Science and Innovation (EWI) (http://www.ewi-vlaanderen.be/en).
#
# http://github.com/hpcugent/easybuild
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
EasyBuild support for BLAT

@author: Andreas Panteli (The Cyprus Institute)
@author: Thekla Loizou (The Cyprus Institute)
"""
import os

from easybuild.easyblocks.generic.makecp import MakeCp
from easybuild.tools.filetools import run_cmd, mkdir

class EB_BLAT(MakeCp):

	def configure_step(self):
		mkdir("bin")

	def build_step(self, verbose=False):
		"""
		Start the actual build
		- typical: make -j X
		"""

		paracmd = ''
		if self.cfg['parallel']:
			paracmd = "-j %s" % self.cfg['parallel']

                bindir=os.path.join(os.getcwd(), "bin")
		
		cmd = "%s make %s BINDIR=%s" % (self.cfg['premakeopts'], paracmd, bindir)

		(out, _) = run_cmd(cmd, log_all=True, simple=False, log_output=verbose)

		return out
