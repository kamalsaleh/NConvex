


#! @Arguments cdd_cone 
#! @Returns a <C>Cone</C> Object
#! @Description  
#! This function takes a cone which is defined using the &GAP; package <C>CddInterface</C> and converts it to a cone in <C>NConvex</C>
DeclareOperation( "Cone",
                  [ IsCddPolyhedron ] );

#! @Arguments C 
#! @Returns a cdd object
#! @Description  
#! Converts the cone to a cdd object. The operations of <C>CddInterface</C> can then be applied
#! on this convex object.
DeclareAttribute( "ExternalCddCone",  IsCone  );

#! @Arguments C 
#! @Returns an normaliz object
#! @Description  
#! Converts the cone to a normaliz object. The operations of NormalizInterface can then be applied
#! on this convex object.
DeclareAttribute( "ExternalNmzCone",  IsCone );

#! @Arguments P
#! @Returns cdd Object
#! @Description  
#! Converts the polyhedron to a cdd object. The operations of CddInterface can then be applied
#! on this convex object.
DeclareAttribute( "ExternalCddPolyhedron",
                   IsPolyhedron );
#! @Arguments P
#! @Returns normaliz Object
#! @Description  
#! Converts the polyhedron to an normaliz object. The operations of NormalizInterface can then be applied
#! on this convex object.
DeclareAttribute( "ExternalNmzPolyhedron",
                   IsPolyhedron );

#! @Arguments P 
#! @Returns a CddPolyhedron
#! @Description  
#! Converts the polytope to a CddPolyhedron. The operations of CddInterface can then be applied
#! on this polyhedron.
DeclareAttribute( "ExternalCddPolytope",
                    IsPolytope );

