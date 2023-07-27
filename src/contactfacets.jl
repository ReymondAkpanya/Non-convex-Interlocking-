#TODO change all "faces" into "facets"
function help_norm(vec::Vector)
    return sqrt(sum(map(i->i^2,vec))) 
end
function help_crossProduct(A::Vector,B::Vector)
    CP=[A[2]*B[3]-A[3]*B[2],
	A[3]*B[1]-A[1]*B[3],
	A[1]*B[2]-A[2]*B[1]];
	return CP;
end

function ComputeContactFacets(poly1::Polyhedron,poly2::Polyhedron)
    res=[]
    coordinates1=poly1.verts
    coordinates2=poly2.verts
    for facet1 in poly1.facets
        for facet2 in poly2.facets
	    intersect=[]
	    print(facet1, facet2,"\n")
	    v=coordinates1[facet1[1]]
	    v1=coordinates1[facet1[2]]-coordinates1[facet1[1]]
	    v2=coordinates1[facet1[3]]-coordinates1[facet1[1]]
	    v3=help_crossProduct(v1,v2)
	    w=coordinates2[facet2[1]]
	    w1=coordinates2[facet2[2]]-coordinates2[facet2[1]]
	    w2=coordinates2[facet2[3]]-coordinates2[facet2[1]]
	    mat=(Matrix([[v1[1], v1[2] ,v1[3]] [v2[1] ,v2[2] ,v2[3]] [v3[1], v3[2] ,v3[3]] ]))^(-1)
	    sol=mat*(w-v)
	    if sol[3]^2<=0.0001
		face1=map(i->coordinates1[i],facet1)
		face2=map(i->coordinates2[i],facet2)
		push!(face1,face1[1])
		push!(face2,face2[1])
		temp=map(i->inpolygon3d(coordinates1[i], face2,1e-8),facet1)
		temp2=map(i->inpolygon3d(coordinates2[i], face1,1e-8),facet2)
		print("temp",temp,temp2,"\n") 		
		
		intersection=[]
		for i in [1,2,3] 
		    if temp[i]^2==1
			push!(intersection,coordinates1[facet1[i]])
		    end
		    if temp2[i]^2==1
			push!(intersection,coordinates2[facet2[i]])
		    end
		end
		print("a",intersection,"\n")
		if map(i->i^2,temp) != [1,1,1] && map(i->i^2,temp2) != [1,1,1]
		    for e1 in EdgesOfFace(poly1,facet1)
		        print("e1=",e1,"\n")
		        v=coordinates1[e1[1]]
			if inpolygon3d(v, face2,1e-8)^2==1 
			    push!(intersection,v) 
			end
			v1=coordinates1[e1[2]]-coordinates1[e1[1]]
			print("v1=",v1,"\n")
			for e2 in EdgesOfFace(poly2,facet2)
			    print("e2=",e2,"\n")
			    w=coordinates2[e2[1]]
			    w1=coordinates2[e2[2]]-coordinates2[e2[1]]
			    n=help_crossProduct(v1,w1)
			    print("n=",n," ",help_norm(n) ,"\n")
			    if help_norm(n)>=1e-8
				mat=(Matrix([[v1[1] ,v1[2], v1[3]] [-w1[1] ,-w1[2], -w1[3]] [n[1] ,n[2], n[3]] ]))^(-1) 
				sol=mat*(w-v)
				print("sol=",sol,"\n")
				if sol[3]^2<1e-8 && sol[1]>=0. && sol[1]<=1. && sol[2]>=0. && sol[2]<=1.
				    push!(intersection,v+sol[1]*v1)
				    print("push\n")
				end
			    end
			end
		    end
		end
		print(intersection,"intersection\n")
		inter=[]
		for i in 1:length(intersection)
		    bool=true	
		    for j in 1:length(inter)   
			if help_norm(intersection[i]-inter[j])<=1e-8
			    bool=false
			end
		    end
		    if bool 
			push!(inter,intersection[i])
		    end
		end
		#if length(inter) >= 3 
		    push!(res,inter)
		#end
	    end
	end
    end
    return res
end


function testComputeContactFacets(num::Int)

    if num==1 
	surf1=Polyhedron([[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
	surf2=Polyhedron([[0.1, 0.1, 0.0], [0.1, 0.2, 0.0], [0.2, 0.1, 0.0]], [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
    elseif num==2 
	surf1=Polyhedron([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]], [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
	surf2=Polyhedron([[0.0, 0.5, 0.0], [1.0, 0.5, 0.0], [0.5, -1.0, 0.0]], [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
    elseif num==3
	surf1=Polyhedron([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.8, 0.5, 0.0]], [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
	surf2=Polyhedron([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.2, 0.5, 0.0]], [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
    elseif num==4 
	surf1=Polyhedron([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]], [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
	surf2=Polyhedron([[1, 1, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]], [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
    end
return ComputeContactFacets(surf1,surf2)

end

