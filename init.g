# SPDX-License-Identifier: GPL-2.0-or-later
# NConvex: A Gap package to perform polyhedral computations
#
# Reading the declaration part of the package.
#

ReadPackage( "NConvex", "gap/ConvexObject.gd" );
ReadPackage( "NConvex", "gap/Fan.gd" );
ReadPackage( "NConvex", "gap/Cone.gd" );
ReadPackage( "NConvex", "gap/Polytope.gd" );
ReadPackage( "NConvex", "gap/Polyhedron.gd" );
ReadPackage( "NConvex", "gap/ZSolve.gd" );

if IsPackageMarkedForLoading( "CddInterface", ">= 2020.06.24" ) and
    IsPackageMarkedForLoading( "NormalizInterface", ">= 1.1.0" ) then
    ReadPackage( "NConvex", "gap/InterfaceToCddNmz.gd" );
fi;
