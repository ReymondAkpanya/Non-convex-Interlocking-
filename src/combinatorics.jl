include("Polyhedron.jl")

#########################################################################
#  TODO: Methode, welche die Facets von Polyeder so umschreibt, dass aufeinanderfolgende Vertizes in Facet durch Kante verbunden sind (siehe bspw. Cuboctahedron, Facet [1,3,5,7])
#########################################################################

#########################################################################
################################## combinatorics
##########################################################################

function NumberOfVertices(poly::Polyhedron)
    return length(poly.verts)
end

function NumberOfEdges(poly::Polyhedron)
    return length(poly.edges)
end

function NumberOfFaces(poly::Polyhedron)
    return length(poly.facets)
end

function EulerCharacteristic(poly::Polyhedron)
    return NumberOfVertices(poly)-NumberOfEdges(poly)+NumberOfFaces(poly)
end


#### EdgesOfFaces

function EdgesOfFaces(poly::Polyhedron)
    res=[]
    return map(i->EdgesOfFace(poly.facets[i]),1:NumberOfFaces(poly))
end

function EdgesOfFace(poly::Polyhedron,face::Vector{Int})
    res=[]
    for i in 1:length(face)-1
	push!(res,[face[i],face[i+1]])
    end
    push!(res,[face[length(face)],face[1]])
    return res
end

function EdgesOfFace(poly::Polyhedron,face::Int)
	return EdgesOfFace(poly,poly.facets[face])
end


####FacesOfEdges
function FacesOfEdges(poly::Polyhedron) 
    return map(i->FacesOfEdge(poly,i),1:NumberOfEdges(poly))
end

function FacesOfEdge(poly::Polyhedron,edge::Vector{<:Int})
    res=[]
    faces=poly.facets
    for face in faces 
        if edge[1] in face && edge[2] in face
	    p1=help_Position(face,edge[1])
	    p2=help_Position(face,edge[2])
	    if p1!= false && p2 != false 
	        if abs(p1-p2)==1 || abs(p1-p2)==length(face)-1
		    push!(res,face)
	        end
	    end
	end
    end
    return res
end

function FacesOfEdge(poly::Polyhedron, edge::Int)
	return FacesOfEdge(poly,poly.edges[edge])
end

#### FacesOfVertices
function FacesOfVertices(poly::Polyhedron )
    res=[]
    for v in 1:NumberOfVertices(poly)
	push!(res,FacesOfVertex(poly,v))
    end
    return res
end

function FacesOfVertex(poly::Polyhedron,v::Int)
    res=[]
    faces=poly.facets
    for f in faces 
        if v in f 
            push!(res,f)
        end
    end
    return res
end

####EdgesOfVertices
function EdgesOfVertex(poly::Polyhedron, v::Int)
    res=[]
    for edge in poly.edges 
	if v in edge
	    push!(res,edge)
	end
    end
    return res
end

function EdgesOfVertices(poly::Polyhedron)
    return map(i->EdgesOfVertex(poly,i),1:NumberOfVertices(poly))
end


####FaceDegreesOfVertices
function FaceDegreesOfVertices(poly::Polyhedron)
    fov=FacesOfVertices(poly)
    return map(i->length(i),fov)
end 

function FaceDegreeOfVertex(poly::Polyhedron,v::Int)
    return length(FacesOfVertex(poly,v))
end

##### Boundaery Edges
function  IsBoundaryEdge(poly::Polyhedron,edge::Vector{<:Int})
    foe=FacesOfEdge(poly,edge)
    return length(foe)==1 
end

function  IsBoundaryEdge(poly::Polyhedron,edge::Int)
    foe=FacesOfEdge(poly,edge)
    return length(foe)==1 
end
function BoundaryEdges(poly::Polyhedron)
    res=[]
    edges=poly.edges
    for edge in edges 
	if IsBoundaryEdge(poly,edge) 
	    push!(res, edge)
	end
    end
    return res
end

####InnerEdges

function  IsInnerEdge(poly::Polyhedron,edge::Int)
    foe=FacesOfEdge(poly,edge)
    return length(foe)==2 
end

function  IsInnerEdge(poly::Polyhedron,edge::Vector{<:Int})
    foe=FacesOfEdge(poly,edge)
    return length(foe)==2 
end

function InnerEdges(poly::Polyhedron)
    res=[]
    edges=poly.edges
    for edge in edges 
	if IsInnerEdge(poly,edge) 
	    push!(res, edge)
	end
    end
    return res
end


#### RamifiedEdges
function IsRamifiedEdge(poly::Polyhedron,edge::Vector{<:Int})
    foe=FacesOfEdge(poly,edge)
    return length(foe)>=3 
end
function IsRamifiedEdge(poly::Polyhedron,edge::Int)
    foe=FacesOfEdge(poly,edge)
    return length(foe)>=3 
end
function RamifiedEdges(poly::Polyhedron)
    res=[]
    edges=poly.edges
    for edge in edges
	if IsRamifiedEdge(poly,edge) 
	    push!(res, edge)
	end
    end
    return res
end

###IsClosedSurface
function IsClosedSurface(poly::Polyhedron)
    return poly.edges==InnerEdges(poly)
end


#### IsCOnnected 
function IsConnected(poly::Polyhedron)
    faces=poly.facets
    visitedFaces=[faces[1]]
    visitedEdges=EdgesOfFace(poly,faces[1])
    faces=filter(f->f!=faces[1],faces)
    tempFaces=0
    while tempFaces != visitedFaces 
	tempF=visitedFaces
	tempE=visitedEdges
	for edge in tempE
	    foe=FacesOfEdge(poly,edge)
	    for f in foe 
		if !(f in visitedFaces)
		    push!(tempF,f)
		    for edge2 in EdgesOfFace(poly,f) 
			if !(edge2 in visitedEdges) && !([edge2[2],edge2[1]] in visitedEdges) 
			    push!(tempE,edge2)
			end
		    end
		end
	    end
	end
	tempFaces=visitedFaces
	visitedFaces=tempF
	visitedEdges=tempE
    end
    return NumberOfFaces(poly)==length(visitedFaces)
end






#### BoundaryVertex
function IsBoundaryVertex(poly::Polyhedron,v::Int)
    eov=EdgesOfVertex(poly,v)
    boundary=[]
    for edge in eov
        if IsBoundaryEdge(poly,edge) 
	    push!(boundary,edge)
	end
    end
    return boundary!=[]
end

function BoundaryVertices(poly::Polyhedron)
    filter(v->IsBoundaryVertex(poly,v),1:NumberOfVertices(poly))
end

####InnerVertex
function IsInnerVertex(poly::Polyhedron,v::Int)
    eov=EdgesOfVertex(poly,v)
    boundary=[]
    for edge in eov
        if IsBoundaryEdge(poly,edge) 
	    push!(boundary,edge)
	end
    end
    return boundary==[]
end

function InnerVertices(poly::Polyhedron)
    filter(v->IsInnerVertex(poly,v),1:NumberOfVertices(poly))
end

## Orientation
function OrientatePolyhedron(poly::Polyhedron)
    faces=poly.facets
    visitedFaces=[faces[1]]
    res=[faces[1]]
    visitedEdges=EdgesOfFace(poly,faces[1])
    faces=filter(f->f!=faces[1],faces)
    tempFaces=0
    while NumberOfFaces(poly) != length(visitedFaces) 
	tempF=visitedFaces
	tempE=visitedEdges
	for edge in tempE
	    foe=FacesOfEdge(poly,edge)
	    for f in foe 
		if !(f in visitedFaces)
		    push!(tempF, f)
		    if edge in EdgesOfFace(poly,f)
			ff=help_OrientateFace(poly,f,edge)
		    else
			ff=f
		    end
		    push!(res,ff)
		    for edge2 in EdgesOfFace(poly,ff) 			
			if !(edge2 in visitedEdges) && !([edge2[2],edge2[1]] in visitedEdges) 
			    push!(tempE,edge2)
			end
		    end
		end
	    end
	end
	visitedFaces=tempF
	visitedEdges=tempE
    end
    return Polyhedron(poly.verts,poly.edges,res)
end

####NeighbourFaceByEdge

function NeighbourFaceByEdge(poly::Polyhedron,face::Vector{<:Int},edge::Vector{<:Int})
   res=[]
   for f in poly.facets
        if edge[1] in f && edge[2] in f 
	    p1=help_Position(face,edge[1])
	    p2=help_Position(face,edge[2])
	    if (abs(p1-p2)==1 || abs(p1-p2)==length(face)-1) && f !=face
		push!(res,f)
	    end
	end 
   end
   return res
end

function NeighbourFaceByEdge(poly::Polyhedron,face::Int,edge::Int)
	return NeighbourFaceByEdge(poly,poly.facets[face],poly.edges[edge])
end

####Disjoint Union
function DisjointUnion(poly::Polyhedron,poly2::Polyhedron)
    num=NumberOfVertices(poly)
    p=CopyPolyhedron(poly)

    for v in poly2.verts
	push!(p.verts,v)
    end
    for edge in poly2.edges
	push!(p.edges,[num+edge[1],num+edge[2]])
    end
    for face in poly2.facets
	push!(p.facets,[num+face[1],num+face[2],num+face[3]])
    end
    return p
end

## face defining vertex cycles
function FaceDefiningVertexCycles(poly::Polyhedron)
  res=[]
  for facet in poly.facets 
    faceCycle=[facet[1]]
    len=1
    while len != length(facet)
      v=faceCycle[len]
      eov=EdgesOfVertex(poly,v)
      temp=filter(v2->[v,v2] in eov || [v2,v] in eov, setdiff(facet,faceCycle))
      push!(faceCycle,temp[1])
      len=len+1
    end    
    push!(res,faceCycle)
  end
  return res
end

################helperfunctions

function help_Position(l::Vector{Vector{<:Int}},x::Vector{<:Int})
    x2=[]
    for i in 1:length(x)
        push!(x2,x[length(x)-i+1])
    end
    for i in 1:length(l)
	if l[i]==x || l[i]==x2
	    return i
	end
    end
    return false
end

function help_Position(l::Vector{Vector{Int}},x::Int)
    for i in 1:length(l)
	if l[i]==x 
	    return i
	end
    end
    return false
end

function help_Position(l::Vector{<:Int},x::Int)
    for i in 1:length(l)
	if l[i]==x 
	    return i
	end
    end
    return false
end

function help_OrientateFace(poly::Polyhedron,face::Vector{<:Int},edge::Vector{<:Int})
    orientFace=[]
    p1=help_Position(face,edge[1])
    p2=help_Position(face,edge[2])
    if p1-p2==-1 || p1-p2==length(face)-1
        for i in 1:length(face)
	    push!(orientFace,face[length(face)-i+1])
	end 
    else
        orientFace=face
    end

    return Vector{Int}(orientFace)

end

function CopyPolyhedron(poly::Polyhedron)
    vertices=map(i->copy(poly.verts[i]),1:NumberOfVertices(poly))
    edges=map(i->copy(poly.edges[i]),1:NumberOfEdges(poly))
    facets=map(i->copy(poly.facets[i]),1:NumberOfFaces(poly))
    return Polyhedron(vertices,edges,facets)

end
