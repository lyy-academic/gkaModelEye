function rayMesh=get_mesh_parallel_rays(rayDirectionVec,centerPosition,radius,numRepeat)

% get_mesh_parallel_rays.m

source_x=centerPosition(1);
source_ys=linspace(-radius,radius,numRepeat)+centerPosition(2);
source_zs=linspace(-radius,radius,numRepeat)+centerPosition(3);



parallel_rays={};
ray_ptr=0;
for i=1:size(source_ys,2)
    for j=1:size(source_zs,2)
        ray_ptr=ray_ptr+1;
        parallel_rays{ray_ptr}=[source_x,source_ys(i),source_zs(j);rayDirectionVec];
    end
end

rayMesh=parallel_rays;


end % get_mesh_parallel_rays