clear;clc;

addpath(genpath(pwd))


sg_1center=createSceneGeometry()
sg_1center.eye.rotationCenters.ele=[-13,0,0];
sg_1center.eye.rotationCenters.azi=[-13,0,0];
sg_1center.eye.rotationCenters.tor=[-13,0,0];

sg_1center.cameraPosition.translation=[0,0,30]';
focal_length=900;
sg_1center.cameraIntrinsic.matrix(1,1)=focal_length;
sg_1center.cameraIntrinsic.matrix(2,2)=focal_length;

% read the stored pose set
poses_set="8directions.csv";
poses=csvread("./poses/"+poses_set,1);


% use linear translation model
translation_method=@bidirectional_linear;
translation_params.direction_azi=[0,1,0];
translation_params.direction_ele=[0,0,1];
translation_params.factor_azi=2;
translation_params.factor_ele=-2;

% use declining translation model
% translation_method=@bidirectional_declining_sin;
% translation_params.direction_azi=[0,1,0];
% translation_params.direction_ele=[0,0,1];
% translation_params.factor_azi=3;
% translation_params.factor_ele=-3;



out_ellipse=[];
rotation_centers_azi=[];
rotation_centers_ele=[];
azi_real=[];
ele_real=[];


for i=1:length(poses)
    azi_real(i)=poses(i,1);
    ele_real(i)=poses(i,2);
    [elli_params_temp,~,image_points,~,~,~,point_labels]=projectModelEye_translation(poses(i,:),sg_1center,'nStopPerimPoints',100,'fullEyeModelFlag',true,'addPseudoTorsion',true,'translation_method',translation_method,'translation_params',translation_params);
    out_ellipse(i,:)=ellipse_transparent2ex(elli_params_temp);
    rotation_centers_azi(i,:)=image_points(find(strcmp(point_labels, 'aziRotationCenter')),:);
    rotation_centers_ele(i,:)=image_points(find(strcmp(point_labels, 'eleRotationCenter')),:);
end

% pure rotation
% cursor=size(out_ellipse,1);
% for i=1:length(poses)
%     [elli_params_temp,~,image_points,~,~,~,point_labels]=projectModelEye(poses{i},sg_1center,'nStopPerimPoints',100,'fullEyeModelFlag',true);
%     out_ellipse(cursor+i,:)=ellipse_transparent2ex(elli_params_temp);
%     rotation_centers_azi(cursor+i,:)=image_points(find(strcmp(point_labels, 'aziRotationCenter')),:);
%     rotation_centers_ele(cursor+i,:)=image_points(find(strcmp(point_labels, 'eleRotationCenter')),:);
% end

names={'cx','cy','a','b','theta','rot_center_azi_x','rot_center_azi_y','rot_center_ele_x','rot_center_ele_y','azi_real','ele_real'};
table2write=table(out_ellipse(:,1),out_ellipse(:,2),out_ellipse(:,3),out_ellipse(:,4),out_ellipse(:,5),rotation_centers_azi(:,1),rotation_centers_azi(:,2),rotation_centers_ele(:,1),rotation_centers_ele(:,2),azi_real',ele_real','VariableNames',names)

write_path="./outputs/";
% writematrix(ellipse_transparent2ex(pupil_ellipse_params),write_path+"test.csv")
writetable(table2write,write_path+poses_set)
