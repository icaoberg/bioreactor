# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

set title "Membrane-coupled bioreactors (MBRs)"
set term png
set output "bioreactor.png"
set xlabel "Time"
set ylabel "Concentration"
set grid
set timestamp
set key top right
set xr [0.0:50.0]
plot "bioreactor.dat" using 1:2 with lines title 'substrate', "bioreactor.dat" using 1:3 with lines title 'biomass'
replot
quit
