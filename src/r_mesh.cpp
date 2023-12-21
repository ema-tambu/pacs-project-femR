// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "include/r_mesh.h"

using cpp_2d_domain = R_Mesh<2,2>;
RCPP_MODULE(cpp_2d_domain) {
    Rcpp::class_<R_Mesh<2,2>>("cpp_2d_domain")
      .constructor<Rcpp::List>()
      .method("nodes"    , &R_Mesh<2,2>::nodes    )
      .method("elements" , &R_Mesh<2,2>::elements )
      .method("neighbors", &R_Mesh<2,2>::neighbors)
      .method("boundary" , &R_Mesh<2,2>::boundary );
}
