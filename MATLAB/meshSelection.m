function [ sampledVertices, sampledFaces, interVertices, interFaces, inter2Faces, vertNumbers] = meshSelection( vertices,faces, boundaries)
%Takes in large set of vertices and faces and returns a bounded portion of them
%   Vertices needs to be 3 columns x,y,z
%   Faces need to be 3 columns of faces based on Vertices input
%   Boundaries need to be [xlower,xupper,ylower,yupper,zlower,zupper]
vertx=vertices(:,1);
verty=vertices(:,2);
vertz=vertices(:,3);
boundx = boundaries(1:2);
boundy = boundaries(3:4);
boundz = boundaries(5:6);
interVertices = zeros(length(vertx),1);
interFaces = zeros(length(faces),1);

%Loop gives 1 or 0 if the vertex is within the boundary
for i=1:length(vertx)
    if (vertx(i)>=boundx(1)) && (vertx(i)<=boundx(2))
        if verty(i) >= boundy(1) && verty(i) <= boundy(2)
            if vertz(i) >= boundz(1) && vertz(i) <= boundz(2)
                interVertices(i)=1;
            else
            end
        else
        end
    else
    end  
end
%Output1
sampledVertices=vertices(find(interVertices),:);
vertNumbers = find(interVertices);

for i=1:length(faces)
    if sum(faces(i,1) == vertNumbers)> 0
        if sum(faces(i,2) == vertNumbers) > 0 
            if sum(faces(i,3) == vertNumbers) > 0
                interFaces(i)=1;
            else
            end
        else
        end
    else
    end
end   

inter2Faces=faces(find(interFaces),:);
for i=1:length(inter2Faces)
    a(1) = find(vertNumbers==inter2Faces(i,1));
    a(2) = find(vertNumbers==inter2Faces(i,2));
    a(3) = find(vertNumbers==inter2Faces(i,3));
    sampledFaces(i,:)=a;
end
%If statements to find which faces use only those vertices


end
