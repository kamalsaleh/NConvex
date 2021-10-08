


##
InstallMethod( RayGenerators,
               [ IsCone ],
               
  function( cone )
#   local nmz_cone, l, r;
  
   return Cdd_GeneratingRays( ExternalCddCone( cone ) );
  
#   nmz_cone := ExternalNmzCone( cone );
  
#   r := NmzGenerators( nmz_cone );
  
#   l := NmzMaximalSubspace( nmz_cone );
  
#   return Concatenation( r, l, -l );

end );

##
InstallMethod( DefiningInequalities,
               [ IsCone ],
               
  function( cone )
  local inequalities, new_inequalities, equalities, i, u; 
  
  inequalities:= ShallowCopy( Cdd_Inequalities( ExternalCddCone( cone ) ) );
  
  equalities:= ShallowCopy( Cdd_Equalities( ExternalCddCone( cone ) ) );
  
  for i in equalities do 
  
       Append( inequalities, [ i,-i ] );
       
  od;
    
  new_inequalities:= [ ];
    
  for i in inequalities do 
  
       u:= ShallowCopy( i );
       
       Remove( u , 1 );
       
       Add(new_inequalities, u );
       
  od;
  
  return new_inequalities; 
    
end );

##
InstallMethod( EqualitiesOfCone,
               "for external Cone",
               [ IsCone ],
               
  function( cone )
  local equalities, new_equalities, u, i;
  
    equalities:= Cdd_Equalities( ExternalCddCone( cone ) );
    
    new_equalities:= [ ];
    
  for i in equalities do 
  
       u:= ShallowCopy( i );
       
       Remove( u , 1 );
       
       Add(new_equalities, u );
       
  od;
  
  return new_equalities;

end );

##
InstallMethod( Dimension, 
               "for cones",
               [ IsCone ],
  function( cone )
 
  return Cdd_Dimension( ExternalCddCone( cone ) );
  
end );

##
InstallMethod( RaysInFacets,
               " for cones",
               [ IsCone ],
               
  function( cone )
  local external_cone, list_of_facets, generating_rays, list, current_cone, current_list, current_ray_generators, i;  
  
    external_cone := Cdd_H_Rep ( ExternalCddCone ( cone ) );
    
    list_of_facets:= Cdd_Facets( external_cone );
    
    generating_rays:= RayGenerators( cone );
    
    list:= [ ];
    
    for i in list_of_facets do
    
      current_cone := Cdd_ExtendLinearity( external_cone, i );
      
      current_ray_generators := Cdd_GeneratingRays( current_cone ) ;
      
      current_list:= List( [1..Length( generating_rays )], 
                           
                           function(j)

                             if generating_rays[j] in Cone( current_cone ) then
                                return 1;
                             else 
                                return 0;
                             fi;
                           
                           end );
                           
      Add( list, current_list );
      
    od;
      
return list;
    
end );

##
InstallMethod( RaysInFaces,
               " for cones",
               [ IsCone ],
               
  function( cone )
  local external_cone, list_of_faces, generating_rays, list, current_cone, current_list, current_ray_generators, i,j;  
  
    external_cone := Cdd_H_Rep( ExternalCddCone ( cone ) );
    
    list_of_faces:= Cdd_Faces( external_cone );
    
    generating_rays:= RayGenerators( cone );
    
    list:= [ ];
    
    for i in list_of_faces do
      
      if i[ 1 ] <> 0 then
        
        current_cone := Cdd_ExtendLinearity( external_cone, i[ 2 ] );
       
        current_ray_generators := Cdd_GeneratingRays( current_cone ) ;
            
        current_list:= List( [ 1 .. Length( generating_rays ) ], 
        
                                function(j)

                                  if generating_rays[j] in Cone( current_cone ) then
                                        return 1;                        
                                  else
                                        return 0;
                                  fi;
                                
                                end );

        Add( list, current_list );

      fi;

   od;

return list;

end );

##
InstallMethod( FVector,
               "for cones",
               [ IsCone ],
  function( cone )
    local external_cone, faces;

    external_cone := Cdd_H_Rep( ExternalCddCone( cone ) );
  
    faces := Cdd_Faces( external_cone );

    return List( [ 1 .. Dimension( cone ) ], 
                i -> Length( PositionsProperty( faces, face -> face[ 1 ] = i ) ) );
end );

##
InstallMethod( IntersectionOfCones,
               "for homalg cones",
               [ IsCone, IsCone ],
               
  function( cone1, cone2 )
    local cone, ext_cone;
    
    if not Rank( ContainingGrid( cone1 ) )= Rank( ContainingGrid( cone2 ) ) then
        
        Error( "cones are not from the same grid" );
        
    fi;
    
    ext_cone := Cdd_Intersection( ExternalCddCone( cone1), ExternalCddCone( cone2 ) );
    
    cone := Cone( ext_cone );
    
    SetContainingGrid( cone, ContainingGrid( cone1 ) );
    
    return cone;
    
end );


#######################################
##
## Methods to construct external cones
##
#######################################

##
InstallMethod( ExternalCddCone, 
               [ IsCone ], 
               
   function( cone )
   
   local list, new_list, number_of_equalities, linearity, i, u ;
   
   new_list:= [ ];
   if IsBound( cone!.input_rays ) and Length( cone!.input_rays )= 1 and IsZero( cone!.input_rays ) then
   
      new_list:= [ Concatenation( [ 1 ], cone!.input_rays[ 1 ] ) ];
      
      return Cdd_PolyhedronByGenerators( new_list );
      
   fi;
   
   if IsBound( cone!.input_rays ) then 
   
      list := cone!.input_rays;
      
      for i in [1..Length( list ) ] do 
          
          u:= ShallowCopy( list[ i ] );
          
          Add( u, 0, 1 );
          
          Add( new_list, u );
      
      od;
      
      return Cdd_PolyhedronByGenerators( new_list );
   
   fi;
   
   
   if IsBound( cone!.input_equalities ) then
   
      list := StructuralCopy( cone!.input_equalities );
      
      number_of_equalities:= Length( list );
      
      linearity := [1..number_of_equalities];
      
      Append( list, StructuralCopy( cone!.input_inequalities ) );
      
      for i in [1..Length( list ) ] do 
      
          u:= ShallowCopy( list[ i ] );
          
          Add( u, 0, 1 );
          
          Add( new_list, u );
      
      od;
      
      return Cdd_PolyhedronByInequalities( new_list, linearity );
   
   else 
   
      list:= StructuralCopy( cone!.input_inequalities );
      
      for i in [1..Length( list ) ] do 
          
          u:= ShallowCopy( list[ i ] );
          
          Add( u, 0, 1 );
          
          Add( new_list, u );
          
      od;
      
      return Cdd_PolyhedronByInequalities( new_list );
   
   fi;
   
end );

##
InstallMethod( ExternalNmzCone, 
              [ IsCone ],
  function( cone )
  local a, list, equalities, i;
  
  list:= [];
   
   if IsBound( cone!.input_rays ) then 
   
        list := StructuralCopy( cone!.input_rays );
        
        if ForAll( Concatenation( list ), IsInt ) then

            return ValueGlobal( "NmzCone" )( [ "integral_closure", list ] );
        
        else

            a := DuplicateFreeList( List( Concatenation( list ), l -> DenominatorRat( l ) ) );
            
            list := Lcm( a ) * list;

            return ValueGlobal( "NmzCone" )( [ "integral_closure", list ] );

        fi;

   fi;
   
   list:= StructuralCopy( cone!.input_inequalities );
   
   if IsBound( cone!.input_equalities ) then
      
      equalities:= StructuralCopy( cone!.input_equalities );
      
      for i in equalities do
      
          Append( list, [ i, -i ] );
          
      od;
      
   fi;
      
    if ForAll( Concatenation( list ), IsInt ) then

        return ValueGlobal( "NmzCone" )( ["inequalities", list ] );
        
    else

        a := DuplicateFreeList( List( Concatenation( list ), l -> DenominatorRat( l ) ) );
            
        list := Lcm( a ) * list;

        return ValueGlobal( "NmzCone" )( ["inequalities", list ] );

    fi;

    Error( "The cone should be defined by vertices or inequalities!" );
    
end );

##
InstallMethod( Cone, 
              "Construct cone from Cdd cone",
              [ IsCddPolyhedron ],
              
   function( cdd_cone )
   local inequalities, equalities, 
         new_inequalities, new_equalities, u, i;
   
   if cdd_cone!.rep_type = "H-rep" then 
       
           inequalities:= Cdd_Inequalities( cdd_cone );
           
           new_inequalities:= [ ];
           
           for i in inequalities do 
                
                 u:= ShallowCopy( i );
                
                 Remove( u , 1 );
                
                 Add(new_inequalities, u );
               
           od;
           
           if cdd_cone!.linearity <> [] then
               
                 equalities:= Cdd_Equalities( cdd_cone );
                 
                 new_equalities:= [ ];
                 
                 for i in equalities do 
                    
                     u:= ShallowCopy( i );
                     
                     Remove( u , 1 );
                     
                     Add(new_equalities, u );
                    
                 od;
                 
                 return ConeByEqualitiesAndInequalities( new_equalities, new_inequalities);
                 
           fi;
           
           return ConeByInequalities( new_inequalities );
           
    else 
    
           return ConeByGenerators( Cdd_GeneratingRays( cdd_cone ) );
           
    fi;
    
end );

##
InstallMethod( IsPointed,
                "for homalg cones.",
                [ IsCone ],
                
   function( cone )
     
     return Cdd_IsPointed( ExternalCddCone ( cone ) );
     
end );

##
InstallMethod( InteriorPoint,
                [ IsConvexObject and IsCone ],
    function( cone )
    local point, denominators;
    point := Cdd_InteriorPoint( ExternalCddCone( cone ) );
    denominators := List( point, DenominatorRat );
    if DuplicateFreeList( denominators ) = [ 1 ] then
        return point;
    else
        return Lcm( denominators )*point;
    fi;
end );

##
InstallMethod( HilbertBasis,
               "for cones",
               [ IsCone ],
               
  function( cone )
    local ineq, const;

    if not IsPointed( cone ) then
        
        Error( "Hilbert basis for not-pointed cones is not yet implemented, you can use the command 'LatticePointsGenerators' " );
        
    fi;


    if IsPackageMarkedForLoading( "NormalizInterface", ">=1.1.0" ) then

      return Set( ValueGlobal( "NmzHilbertBasis" )( ExternalNmzCone( cone ) ) );
    
    elif IsPackageMarkedForLoading( "4ti2Interface", ">=2018.07.06" ) then
      
      ineq := DefiningInequalities( cone );

      const := ListWithIdenticalEntries( Length( ineq ), 0 );

      return Set( ValueGlobal( "4ti2Interface_zsolve_equalities_and_inequalities" )( [  ], [  ], ineq, const )[ 2 ]: precision := "gmp" );

    else

      Error( "4ti2Interface or NormalizInterface should be loaded!" );

    fi;
  
end );


##############
#
# Polytopes
#
#############

##
InstallMethod( IsEmpty,
               "for polytopes",
               [ IsPolytope ],
               
  function( polytope )
    
    if IsBound( polytope!.input_points ) and Length( polytope!.input_points ) > 0 then
        
        return false;
        
    elif IsBound( polytope!.input_points ) and Length( polytope!.input_points ) = 0 then
    
        return true;
    
    else 
    
       return Cdd_IsEmpty( ExternalCddPolytope( polytope ) );
       
    fi;
    
end );

##
InstallMethod( InteriorPoint,
                [ IsConvexObject and IsPolytope ],
    function( poly )
    return Cdd_InteriorPoint( ExternalCddPolytope( poly ) );
end );

##
InstallMethod( IsBounded,
               " for external polytopes.",
               [ IsPolytope ],
               
  function( polytope )

  return Length( Cdd_GeneratingRays( ExternalCddPolytope( polytope ) ) ) = 0;
  
end );

##
InstallMethod( ExternalCddPolytope, 
               "for polytopes", 
               [ IsPolytope ],
               
   function( polyt )
   local old_pointlist, new_pointlist, ineqs, i,j;
   
   if IsBound( polyt!.input_points ) and IsBound( polyt!.input_ineqs ) then
        
        Error( "points and inequalities at the same time are not supported\n" );
        
   fi;
    
   if IsBound( polyt!.input_points ) then 
   
       old_pointlist := polyt!.input_points;
       
       new_pointlist:= [ ];
       
       for i in old_pointlist do 
           
           j:= ShallowCopy( i );
           
           Add( j, 1, 1 );
           
           Add( new_pointlist, j );
           
       od;
           
       return Cdd_PolyhedronByGenerators( new_pointlist );
       
   elif  IsBound( polyt!.input_ineqs ) then
    
      ineqs := ShallowCopy( polyt!.input_ineqs );
      
      return Cdd_PolyhedronByInequalities( ineqs );
      
   else 
   
       Error( "something went wrong\n" );
       
   fi;
   
end );

##
InstallMethod( Dimension, 
               "for polytopes",
               [ IsPolytope ],
  function( polytope )

    return Cdd_Dimension( ExternalCddPolytope( polytope ) );
    
end );

##
InstallMethod( VerticesOfPolytope,
               "for polytopes",
               [ IsPolytope ],
               
  function( polyt )
    
    return Cdd_GeneratingVertices( ExternalCddPolytope( polyt ) );
    
end );

##
InstallMethod( FacetInequalities,
               " for external polytopes",
               [ IsExternalPolytopeRep ],
               
  function( polyt )
    
    return Cdd_Inequalities( ExternalCddPolytope( polyt ) );
    
end );

##
InstallMethod( EqualitiesOfPolytope,
               "for external polytopes",
               [ IsPolytope ],
               
  function( polyt )
    
    return Cdd_Equalities( ExternalCddPolytope( polyt ) );
    
end );

##
InstallMethod( FVector,
        "for polytopes",
        [ IsPolytope ],
    function( polyt )
      local external_polytope, faces;
      
      external_polytope := Cdd_H_Rep( ExternalCddPolytope( polyt ) );
      
      faces := Cdd_Faces( external_polytope );
      
      return List( [ 0 .. Dimension( polyt ) - 1 ], 
                i -> Length( PositionsProperty( faces, face -> face[ 1 ] = i ) ) );

end );

##
InstallMethod( IntersectionOfPolytopes,
               "for homalg cones",
               [ IsPolytope, IsPolytope ],
               
  function( polyt1, polyt2 )
    local polyt, ext_polytope;
    
    if not Rank( ContainingGrid( polyt1 ) ) = Rank( ContainingGrid( polyt2 ) ) then
        
        Error( "polytopes are not of the same dimension" );
        
    fi;
    
    ext_polytope:= Cdd_Intersection( ExternalCddPolytope( polyt1), ExternalCddPolytope( polyt2) ); 
    
    polyt := Polytope( Cdd_GeneratingVertices( ext_polytope) );
    
    SetExternalCddPolytope( polyt, ext_polytope );
    
    SetContainingGrid( polyt, ContainingGrid( polyt1 ) );
    
    SetAmbientSpaceDimension( polyt, AmbientSpaceDimension( polyt1 ) );
    
    return polyt;
    
end );

##
InstallMethod( IsEmpty,
               "for polytopes",
               [ IsPolytope ],
               
  function( polytope )
    
    if IsBound( polytope!.input_points ) and Length( polytope!.input_points ) > 0 then
        
        return false;
        
    elif IsBound( polytope!.input_points ) and Length( polytope!.input_points ) = 0 then
    
        return true;
    
    else 
    
       return Cdd_IsEmpty( ExternalCddPolytope( polytope ) );
       
    fi;
    
end );

##
InstallMethod( InteriorPoint,
                [ IsConvexObject and IsPolytope ],
    function( poly )
    return Cdd_InteriorPoint( ExternalCddPolytope( poly ) );
end );

##
InstallMethod( IsBounded,
               " for external polytopes.",
               [ IsPolytope ],
               
  function( polytope )

  return Length( Cdd_GeneratingRays( ExternalCddPolytope( polytope ) ) ) = 0;
  
end );

##
InstallMethod( ExternalCddPolytope, 
               "for polytopes", 
               [ IsPolytope ],
               
   function( polyt )
   local old_pointlist, new_pointlist, ineqs, i,j;
   
   if IsBound( polyt!.input_points ) and IsBound( polyt!.input_ineqs ) then
        
        Error( "points and inequalities at the same time are not supported\n" );
        
   fi;
    
   if IsBound( polyt!.input_points ) then 
   
       old_pointlist := polyt!.input_points;
       
       new_pointlist:= [ ];
       
       for i in old_pointlist do 
           
           j:= ShallowCopy( i );
           
           Add( j, 1, 1 );
           
           Add( new_pointlist, j );
           
       od;
           
       return Cdd_PolyhedronByGenerators( new_pointlist );
       
   elif  IsBound( polyt!.input_ineqs ) then
    
      ineqs := ShallowCopy( polyt!.input_ineqs );
      
      return Cdd_PolyhedronByInequalities( ineqs );
      
   else 
   
       Error( "something went wrong\n" );
       
   fi;
   
end );

##
InstallMethod( Dimension, 
               "for polytopes",
               [ IsPolytope ],
  function( polytope )
    
    return Cdd_Dimension( ExternalCddPolytope( polytope ) );
    
end );

##
InstallMethod( VerticesOfPolytope,
          "for polytopes",
          [ IsPolytope ],           
  function( polyt )
    
    return Cdd_GeneratingVertices( ExternalCddPolytope( polyt ) );
    
end );

##
InstallMethod( FacetInequalities,
               " for external polytopes",
               [ IsExternalPolytopeRep ],
               
  function( polyt )
    
    return Cdd_Inequalities( ExternalCddPolytope( polyt ) );
    
end );

##
InstallMethod( EqualitiesOfPolytope,
               "for external polytopes",
               [ IsPolytope ],
               
  function( polyt )
    
    return Cdd_Equalities( ExternalCddPolytope( polyt ) );
    
end );

##
InstallMethod( FVector,
        "for polytopes",
        [ IsPolytope ],
    function( polyt )
      local external_polytope, faces;
      
      external_polytope := Cdd_H_Rep( ExternalCddPolytope( polyt ) );
      
      faces := Cdd_Faces( external_polytope );
      
      return List( [ 0 .. Dimension( polyt ) - 1 ], 
                i -> Length( PositionsProperty( faces, face -> face[ 1 ] = i ) ) );

end );

##
InstallMethod( IntersectionOfPolytopes,
               "for homalg cones",
               [ IsPolytope, IsPolytope ],
               
  function( polyt1, polyt2 )
    local polyt, ext_polytope;
    
    if not Rank( ContainingGrid( polyt1 ) ) = Rank( ContainingGrid( polyt2 ) ) then
        
        Error( "polytopes are not of the same dimension" );
        
    fi;
    
    ext_polytope:= Cdd_Intersection( ExternalCddPolytope( polyt1), ExternalCddPolytope( polyt2) ); 
    
    polyt := Polytope( Cdd_GeneratingVertices( ext_polytope) );
    
    SetExternalCddPolytope( polyt, ext_polytope );
    
    SetContainingGrid( polyt, ContainingGrid( polyt1 ) );
    
    SetAmbientSpaceDimension( polyt, AmbientSpaceDimension( polyt1 ) );
    
    return polyt;
    
end );



##############
#
# Polyhedron
#
##############

##
InstallMethod( ExternalCddPolyhedron,
               "for polyhedrons",
               [ IsPolyhedron and HasMainRatPolytope and HasTailCone ],
               
  function( polyhedron )
    local verts, rays;
    
    verts := Vertices( MainRatPolytope( polyhedron ) );
    
    verts := List( verts, i -> Concatenation( [ 1 ], i ) );
    
    rays := RayGenerators( TailCone( polyhedron ) );
    
    rays := List( rays, i -> Concatenation( [ 0 ], i ) );
    
    polyhedron := Concatenation( rays, verts );
    
    polyhedron := Cdd_PolyhedronByGenerators( polyhedron );
    
    return polyhedron;
    
end );

##
InstallMethod( InteriorPoint,
                [ IsConvexObject and IsPolyhedron ],
    function( poly )
    return Cdd_InteriorPoint( ExternalCddPolyhedron( poly ) );
end );

##
InstallMethod( ExternalCddPolyhedron,
               "for polyhedrons with inequalities",
               [ IsPolyhedron ],
               
  function( polyhedron )
    
    if IsBound( polyhedron!.inequalities ) then
        
        if IsEmpty( polyhedron!.inequalities ) then
            
            polyhedron!.inequalities := [ [ 0 ] ];
            
        fi;
        
        return Cdd_PolyhedronByInequalities( polyhedron!.inequalities );
        
    fi;
    
    TryNextMethod();
    
end );

##
InstallMethod( ExternalNmzPolyhedron, 
               [ IsPolyhedron ], 
  function( poly )
    local ineq, new_ineq;
    
    if IsBound( poly!.inequalities ) then 
      
      ineq := poly!.inequalities;
    
    fi;
    
    ineq := DefiningInequalities( poly );
    
    new_ineq :=
      List( ineq,
        function( i )
          local j;
          j:= ShallowCopy( i );
          Add( j, j[ 1 ] );
          Remove(j ,1 );
          return j;
      end );
    
    return ValueGlobal( "NmzCone" )( [ "inhom_inequalities", new_ineq ] );
    
end );

##
InstallMethod( Dimension,
               [ IsPolyhedron ],
   function( polyhedron )
   
   return Cdd_Dimension( ExternalCddPolyhedron( polyhedron ) );
   
end );

##
InstallMethod( VerticesOfMainRatPolytope,
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    local v;
    
    if IsBound( polyhedron!.inequalities ) then
        
        v:= Cdd_GeneratingVertices( ExternalCddPolyhedron( polyhedron ) );
        
        if Length( v ) > 0 then
          
          return v;
          
        else
          
          return [ ListWithIdenticalEntries(AmbientSpaceDimension( polyhedron ), 0 ) ];
          
        fi;
        
    else
      
      return Vertices( MainRatPolytope( polyhedron ) );
      
    fi;
    
end );

##
InstallMethod( RayGeneratorsOfTailCone,
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    
    if IsBound( polyhedron!.inequalities ) then
      
      return Cdd_GeneratingRays( ExternalCddPolyhedron( polyhedron ) );
      
    else
      
      return RayGenerators( TailCone( polyhedron ) );
      
    fi;
    
end );

##
InstallMethod( FVector,
            [ IsPolyhedron ],
    function( polyhedron )
      local external_polyhedron, faces;
      
      external_polyhedron := Cdd_H_Rep( ExternalCddPolyhedron( polyhedron ) );
      
      faces := Cdd_Faces( external_polyhedron );
      
      return List( [ 0 .. Dimension( polyhedron ) - 1 ],
                i -> Length( PositionsProperty( faces, face -> face[ 1 ] = i ) ) );

end );

## Solving linear programs
##
BindGlobal( "SOLVE_LINEAR_PROGRAM_USING_CDD",
  function( P, max_or_min, target_func, constructor )
    local ext_cdd_poly, cdd_linear_program; 
    
    ext_cdd_poly := constructor( P );
    
    cdd_linear_program := Cdd_LinearProgram( ext_cdd_poly, max_or_min, target_func );
    
    return Cdd_SolveLinearProgram( cdd_linear_program );

end );

##
InstallMethod( SolveLinearProgram,
  [ IsPolyhedron, IsString, IsList ],
  function( P, max_or_min, target_func )
    
    return SOLVE_LINEAR_PROGRAM_USING_CDD( P, max_or_min, target_func, ExternalCddPolyhedron );

end );

##
InstallMethod( SolveLinearProgram,
  [ IsPolytope, IsString, IsList ],
  function( P, max_or_min, target_func )
    
    return SOLVE_LINEAR_PROGRAM_USING_CDD( P, max_or_min, target_func, ExternalCddPolytope );

end );


##
InstallMethod( ExternalCddPolyhedron,
               "for polyhedrons",
               [ IsPolyhedron and HasHomogeneousPointsOfPolyhedron ],
               
  function( polyhedron )
    
    return Cdd_PolyhedronByGenerators( HomogeneousPointsOfPolyhedron( polyhedron ) );
    
end );

##
InstallMethod( DefiningInequalities, 
               " for polyhedrons",
               [ IsPolyhedron ], 
               
   function( polyhedron )
   local ineq, eq, ex, d;
   
   ex:= ExternalCddPolyhedron( polyhedron );
   
   ineq := Cdd_Inequalities( ex );
   eq   := Cdd_Equalities( ex );
   
   d:= Concatenation( ineq, eq, -eq );
   
   return d;
 
end );

##
InstallMethod( LatticePointsGenerators,
                 [ IsPolyhedron ],
                 
  function( p )
    local external_poly,nmz_points_in_main_polytope, points_in_main_polytope,
          nmz_hilbert_basis, hilbert_basis, nmz_lineality, lineality, ineq, const;
    
    if IsPackageMarkedForLoading( "NormalizInterface", ">=1.1.0" ) then

      external_poly:= ExternalNmzPolyhedron( p );

      nmz_points_in_main_polytope:= ValueGlobal( "NmzModuleGenerators" )( external_poly );

      points_in_main_polytope:=
        List( nmz_points_in_main_polytope ,
          function( i )
            local j;
            
            j := ShallowCopy( i );
            
            Remove( j, Length( i ) );
            
            return j;
            
          end );
          
      nmz_hilbert_basis:= ValueGlobal( "NmzHilbertBasis" )( external_poly );
      
      hilbert_basis :=
        List( nmz_hilbert_basis ,
          function( i )
            local j;
            
            j := ShallowCopy( i );
            
            Remove( j, Length( i ) );
            
            return j;
            
          end );
          
      nmz_lineality := ValueGlobal( "NmzMaximalSubspace" )( external_poly );
      
      lineality:= List( nmz_lineality,
        function( i )
          local j;
          
          j := ShallowCopy( i );
          
          Remove( j, Length( i ) );
          
          return j;
          
        end );
        
      return [ Set( points_in_main_polytope ), Set( hilbert_basis ), Set( lineality ) ];
      
    elif IsPackageMarkedForLoading( "4ti2Interface", ">=2018.07.06" ) then
     
      ineq := TransposedMat( DefiningInequalities( p ) );
      
      const := -ineq[ 1 ];
      
      ineq := TransposedMat( ineq{ [ 2 .. Length( ineq ) ] } );
      
      return List( ValueGlobal( "4ti2Interface_zsolve_equalities_and_inequalities" )( [  ], [  ], ineq, const : precision := "gmp" ), Set );
      
    else
      
      Error( "4ti2Interface or NormalizInterface should be loaded!" );
      
    fi;
    
end );

##
InstallMethod( IsBounded,
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    
    return Length( Cdd_GeneratingRays( ExternalCddPolyhedron( polyhedron ) ) ) = 0;
    
end );

##
InstallMethod( IsNotEmpty,
               "for polyhedrons",
               [ IsPolyhedron ],
               
  function( polyhedron )
    
    return not Cdd_IsEmpty( ExternalCddPolyhedron( polyhedron ) );
    
end );
