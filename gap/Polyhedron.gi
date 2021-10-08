# SPDX-License-Identifier: GPL-2.0-or-later
# NConvex: A Gap package to perform polyhedral computations
#
# Implementations
#



DeclareRepresentation( "IsConvexPolyhedronRep",
                       IsPolyhedron and IsExternalConvexObjectRep,
                       [ ]
                      );

####################################
##
## Types and Families
##
####################################


BindGlobal( "TheFamilyOfPolyhedrons",
        NewFamily( "TheFamilyOfPolyhedrons" , IsPolyhedron ) );

BindGlobal( "TheTypeConvexPolyhedron",
        NewType( TheFamilyOfPolyhedrons,
                 IsPolyhedron and IsConvexPolyhedronRep ) );
                                 
#####################################
##
## Structural Elements
##
#####################################

##
InstallMethod( ContainingGrid,
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    
    if HasTailCone( polyhedron ) then
        
        return ContainingGrid( TailCone( polyhedron ) );
        
    elif HasMainRatPolytope( polyhedron ) then
        
        return ContainingGrid( MainRatPolytope( polyhedron ) );
        
    fi;
    
end );

##
InstallMethod( MainRatPolytope, 
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    
    return Polytope( VerticesOfMainRatPolytope( polyhedron ) );
    
end );

##
InstallMethod( MainPolytope,
              [ IsPolyhedron ],
  
  function( polyhedron )
    local V;

    V := VerticesOfMainRatPolytope( polyhedron );

    if ForAll( V, v -> ForAll( v, IsInt ) ) then

      return MainRatPolytope( polyhedron );

    fi;
  
    return Polytope( LatticePointsGenerators( polyhedron )[ 1 ] );
    
end );

##
InstallMethod( VerticesOfMainPolytope,
              [ IsPolyhedron ], 
              
  function( polyhedron )
    local V;
    
    V := VerticesOfMainRatPolytope( polyhedron );
    
    if ForAll( V, v -> ForAll( v, IsInt ) ) then

      return V;

    fi;

    return Vertices( MainPolytope( polyhedron ) );

end );

##
InstallMethod( TailCone,
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    
  if RayGeneratorsOfTailCone( polyhedron ) <> [ ] then 
  
         return Cone( RayGeneratorsOfTailCone( polyhedron ) );
         
  else 
  
  
         return Cone( [ ListWithIdenticalEntries( AmbientSpaceDimension( polyhedron ), 0 ) ] );
   
  fi;
  
end );

##
InstallMethod( HomogeneousPointsOfPolyhedron,
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    local verts, rays;
    
    verts := Vertices( MainRatPolytope( polyhedron ) );
    
    verts := List( verts, i -> Concatenation( [ 1 ], i ) );
    
    rays := RayGenerators( TailCone( polyhedron ) );
    
    rays := List( rays, i -> Concatenation( [ 0 ], i ) );
    
    polyhedron := Concatenation( rays, verts );
    
    return polyhedron;
    
end );
    
##
InstallMethod( BasisOfLinealitySpace,
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    
    return LinealitySpaceGenerators( TailCone( polyhedron ) );
    
end );

#####################################
##
##  Properties
##
#####################################

##
InstallMethod( IsPointed,
               [ IsPolyhedron ],
               
    function( polyhedron )
      
      return IsPointed( TailCone( polyhedron ) );
      
end );

#####################################
##
## Constructors
##
#####################################

##
InstallMethod( PolyhedronByInequalities,
               "for list of inequalities",
               [ IsList ],
               
  function( inequalities )
    local polyhedron;
    
    polyhedron := rec();
    
    ObjectifyWithAttributes( polyhedron, TheTypeConvexPolyhedron );
    
    polyhedron!.inequalities := inequalities;
    
    if not IsEmpty( inequalities ) then
        
        SetAmbientSpaceDimension( polyhedron, Length( inequalities[ 1 ] ) - 1 );
        
    else
        
        SetAmbientSpaceDimension( polyhedron, 0 );
    fi;
    
    return polyhedron;
    
end );

##
InstallMethod( Polyhedron,
               "for a polytope and a cone",
               [ IsPolytope, IsCone ],
               
  function( polytope, cone )
    local polyhedron;
    
    if not Rank( ContainingGrid( polytope ) )= Rank( ContainingGrid( cone ) ) then
        
        Error( "Two objects are not comparable" );
        
    fi;
    
    polyhedron := rec();
    
    ObjectifyWithAttributes( polyhedron, TheTypeConvexPolyhedron,
                                          MainRatPolytope, polytope,
                                          TailCone, cone,
                                          ContainingGrid, ContainingGrid( polytope ),
                                          AmbientSpaceDimension, AmbientSpaceDimension( polytope )
                                        );
    
    return polyhedron;
    
end );

##
InstallMethod( Polyhedron,
               "for a polytope and a list",
               [ IsPolytope, IsList ],
               
  function( polytope, cone )
    local polyhedron;
    
    if Length( cone ) > 0 and Length( cone[ 1 ] ) <> AmbientSpaceDimension( polytope ) then
        
        Error( "the two objects are not comparable" );
        
    fi;
    
    polyhedron := rec( );
    
    if Length( cone ) = 0 then
        
        cone := [ List( [ 1 .. AmbientSpaceDimension( polytope ) ], i -> 0 ) ];
        
    fi;
    
    ObjectifyWithAttributes( polyhedron, TheTypeConvexPolyhedron,
                                          MainRatPolytope, polytope,
                                          TailCone, Cone( cone ),
                                          ContainingGrid, ContainingGrid( polytope ),
                                          AmbientSpaceDimension, AmbientSpaceDimension( polytope )
                                        );
    
    return polyhedron;
    
end );


##
InstallMethod( Polyhedron,
               "for a polytope and a cone",
               [ IsList, IsCone ],
               
  function( polytope, cone )
    local polyhedron;
    
    if Length( polytope ) > 0 and Length( polytope[ 1 ] ) <> AmbientSpaceDimension( cone ) then
        
        Error( "the two objects are not comparable" );
        
    fi;
    
    polytope := Polytope( polytope );
    
    SetContainingGrid( polytope, ContainingGrid( cone ) );
    
    polyhedron := rec( );
    
    ObjectifyWithAttributes( polyhedron, TheTypeConvexPolyhedron,
                                          MainRatPolytope, polytope,
                                          TailCone, cone,
                                          ContainingGrid, ContainingGrid( cone ),
                                          AmbientSpaceDimension, AmbientSpaceDimension( cone )
                                        );
    
    return polyhedron;
    
end );

##
InstallMethod( Polyhedron,
               "for a polytope and a cone",
               [ IsList, IsList ],
               
  function( polytope, cone )
    local polyhedron;
    
    if Length( polytope ) > 0 and Length( cone ) > 0 and Length( cone[ 1 ] ) <> Length( polytope[ 1 ] ) then
        
        Error( "two objects are not comparable\n" );
        
    fi;
    
    if Length( polytope ) = 0 then
        
        Error( "The polytope of a polyhedron should at least contain one point!" );
        
    fi;
    
    if Length( cone ) = 0 then
        
        cone := [ List( [ 1 .. Length( polytope[ 1 ] ) ], i -> 0 ) ];
        
    fi;
    
    polyhedron := rec();
    
    ObjectifyWithAttributes( polyhedron, TheTypeConvexPolyhedron,
                                          MainRatPolytope, Polytope( polytope ),
                                          TailCone, Cone( cone ),
                                          AmbientSpaceDimension, Length( polytope[ 1 ] ) 
                                        );
    
    SetContainingGrid( TailCone( polyhedron ), ContainingGrid( MainRatPolytope( polyhedron ) ) );
    
    SetContainingGrid( polyhedron, ContainingGrid( MainRatPolytope( polyhedron ) ) );
    
    return polyhedron;
    
end );

##############################
##
## View & Display
##
##############################

##
InstallMethod( ViewObj,
               "for homalg polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    local str;
    
    Print( "<A" );
    
    if HasIsNotEmpty( polyhedron ) then
        
        if IsNotEmpty( polyhedron ) then
            
            Print( " not empty" );
            
        fi;
    
    fi;
    
    Print( " polyhedron in |R^" );
    
    Print( String( AmbientSpaceDimension( polyhedron ) ) );
    
    if HasDimension( polyhedron ) then
        
        Print( " of dimension ", String( Dimension( polyhedron ) ) );
        
    fi;
    
    Print( ">" );
    
end );

##
InstallMethod( Display,
               "for homalg polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    local str;
    
    Print( "A" );
    
    if HasIsNotEmpty( polyhedron ) then
        
        if IsNotEmpty( polyhedron ) then
            
            Print( " not empty" );
            
        fi;
    
    fi;
    
    Print( " polyhedron in |R^" );
    
    Print( String( AmbientSpaceDimension( polyhedron ) ) );
    
    if HasDimension( polyhedron ) then
        
        Print( " of dimension ", String( Dimension( polyhedron ) ) );
        
    fi;
    
    Print( ".\n" );
    
end );
