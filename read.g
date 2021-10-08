# SPDX-License-Identifier: GPL-2.0-or-later
# NConvex: A Gap package to perform polyhedral computations
#
# Reading the implementation part of the package.
#

ReadPackage( "NConvex", "gap/ConvexObject.gi" );
ReadPackage( "NConvex", "gap/Fan.gi" );
ReadPackage( "NConvex", "gap/Cone.gi");
ReadPackage( "NConvex", "gap/Polytope.gi");
ReadPackage( "NConvex", "gap/Polyhedron.gi");
ReadPackage( "NConvex", "gap/ZSolve.gi");

if IsPackageMarkedForLoading( "CddInterface", ">= 2020.06.24" ) and
    IsPackageMarkedForLoading( "NormalizInterface", ">= 1.1.0" ) then
    ReadPackage( "NConvex", "gap/InterfaceToCddNmz.gi" );
fi;
