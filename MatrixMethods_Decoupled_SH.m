%% Matrix Methods Decoupled SH
% Update 2026/05/23
% written by: Armin Huber, armin.huber@dlr.de
% -------------------------------------------------------------------------
% Calculate SH wave dispersion diagrams for free anisotropic multilayered 
% plates in decoupled cases by using the transfer matrix method (TMM) or 
% the stiffness matrix method (SMM). The dispersion  equation amplitude is 
% scaled to a grayscale. Modal solutions are represented by minima, i.e., 
% by dark shading. The numerical instability of TMM is no problem here 
% because it occurs only below the bulk shear velocity. Because all guided 
% SH modes are above, they are not affected. SMM is stable, but it does not 
% find all modes. Therefore, for computing SH modes, using TMM is the 
% better choice.

% NOTE: This code is valid for wave propagation along axes of symmetry
% only. The layup may only contain single or cross plies and wave 
% propagation must be along or normal to an axis of symmetry. The equations 
% decouple into motion in the sagittal plane x1-x3 (pure Lamb modes) and 
% motion in the shear horizontal direction x2 (pure SH modes). This code 
% computes the SH modes.

% References
% [1] A. H. Nayfeh, "The general problem of elastic wave propagation in
%     multilayered anisotropic media," J. Acoust. Soc. Am. 89(4), 1521–1531
%     (1991).
% [2] A. H. Nayfeh, "The propagation of horizontally polarized shear waves 
%     in multilayered anisotropic media," J. Acoust. Soc. Am. 86(5), 2007–
%     2012 (1989).
% [3] L. Wang and S. I. Rokhlin, "Stable reformulation of transfer matrix 
%     method for wave propagation in layered anisotropic media," 
%     Ultrasonics 39, 413-424 (2001).
% [4] S. I. Rokhlin and L. Wang, "Stable recursive algorithm for elastic
%     wave propagation in layered anisotropic media: Stiffness matrix
%     method," J. Acoust. Soc. Am. 112(3), 822-834 (2002).
% [5] A. M. A. Huber, "Classification of solutions for guided waves in
%     fluid-loaded viscoelastic composites with large numbers of layers,"
%     J. Acoust. Soc. Am. 154(2), 1073–1094 (2023).

%% settings
clear

Name = 'T800M913';
Density = 1550; % (kg/m^3)
C(1,1) = 154; C(1,2) = 3.7; C(1,3) = 3.7; % stiffness (GPa)
              C(2,2) = 9.5; C(2,3) = 5.2;
                            C(3,3) = 9.5;
                                          C(4,4) = 2.15;
                                                         C(5,5) = 4.2;
                                                                       C(6,6) = 4.2;

LayerOrientations = [0 90]; % orientation of the layers in the super layer (deg), e.g., [0 90]
LayerThicknesses = [.125 .125]; % corresponding layer thicknesses (mm), e.g., [.125 .125]

Repetitions = 2; % number of repetitions of the super layer, e.g., 2 for [0 90]2 = [0 90 0 90]
SymmetricSystem = 1; % 1: yes 0: no; if yes, the laminate is mirrored, e.g., [0 90]2s = [0 90 0 90 90 0 90 0]

MatrixMethod = 'TMM'; % select matrix method: 'TMM' or 'SMM' (transfer matrix method or stiffness matrix method)

PropagationAngle = 0; % wave propagation direction with respect to LayerOrientations (deg)
PhaseVelocityLimit = 20; % (m/ms)
FrequencyLimit = 4000; % (kHz)
PhaseVelocitySteps = 500;
FrequencySteps = 500;

%% calculate transformed stiffnesses
C = C*1e9; % stiffness (Pa)
SuperLayerSize = length(LayerOrientations);
Beta = LayerOrientations-PropagationAngle; % angle between the fibers of each layer and the propagation angle
if  any(mod(Beta,90)) % check if the layup supports decoupling and wave propagation is along an axis of symmetry
    errordlg('Unsupported layup or propagation direction. "LayerOrientations" may only contain single layers or cross plies such as 0 and 90 degree layers, and "PropagationAngle" may only be along or normal to a principal axis, e.g., 0 or 90 degrees in this case.','Error');
    return
end
for m = 1:SuperLayerSize % transformed stiffnesses for orthotropic materials for rotation about x3 [1], Appendix A
    s = sind(Beta(m));
    g = cosd(Beta(m));
    c{m}(4,4) = C(4,4)*g^2+C(5,5)*s^2; %#ok<*SAGROW>
    c{m}(6,6) = C(6,6)+(C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2;
end
n = 2.^(0:log2(Repetitions));
Pattern = [0:length(n)-1;n]; % for efficiently using the recursive algorithm on the super layer transfer/stiffness matrices to obtain the global transfer/stiffness matrix
while Pattern(2,end) < Repetitions
    for i = length(n):-1:1
        if  Pattern(2,end)+n(i) <= Repetitions
            Pattern(1:2,end+1) = [i;Pattern(2,end)+n(i)];
            break
        end
    end
end
Pattern = Pattern(1,2:end);
LayerThicknesses = LayerThicknesses/1e3; % (m)
Layers = Repetitions*SuperLayerSize;
Thickness = Repetitions*sum(LayerThicknesses)*1e3;
if  SymmetricSystem %#ok<*UNRCH>
    Layers = 2*Layers; 
    Thickness = 2*Thickness;
end
disp(['Method:    ',MatrixMethod,newline...
      'Layers:    ',num2str(Layers),newline...
      'Thickness: ',num2str(Thickness),' mm'])
warning('off','MATLAB:nearlySingularMatrix') % turn warnings off
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:illConditionedMatrix')

%% calculate dispersion diagram
tic
Frequency = FrequencyLimit/FrequencySteps:FrequencyLimit/FrequencySteps:FrequencyLimit; % generate a range of frequencies
PhaseVelocity = (PhaseVelocityLimit/PhaseVelocitySteps:PhaseVelocityLimit/PhaseVelocitySteps:PhaseVelocityLimit)*1e3; % generate a range of phase velocities
AngularFrequency = 2*pi*Frequency*1e3;
Wavenumber = AngularFrequency./PhaseVelocity';
Size = size(Wavenumber);
Length = numel(Wavenumber);
Y = complex(zeros(Size)); %#ok<*NASGU> % preallocate memory for the dispersion equation amplitude
k2 = reshape(Wavenumber.^2,1,1,Length); % reshape into pages in order to use the pagewise matrix function "pagemtimes"; this enables implicit parallel computing and we avoid slow loops
rw2 = reshape(repmat(Density*AngularFrequency.^2,Size(1),1),1,1,Length);
for m = 1:SuperLayerSize
    k3 = sqrt((rw2-k2*c{m}(6,6))/c{m}(4,4)); % out-of-plane wavenumber component of the pair of bulk waves (SH+-) in [5], Eq. (48)
    D4 = 1i*k3*c{m}(4,4); % shear stress amplitude (sigma23) [5], Eq. (48)
    if  strcmp(MatrixMethod,'TMM')
        G = k3*LayerThicknesses(m);
        CosG = cos(G);
        SinG = sin(G);
        L{m} = [CosG 1i*SinG./D4;1i*SinG.*D4 CosG]; % mth layer transfer matrix [5], Eq. (49)
    elseif strcmp(MatrixMethod,'SMM')
        E = exp(1i*k3*LayerThicknesses(m)); % harmonic term in the out-of-plane direction in [5], Eqs. (46)-(47)
        E2 = E.^2;
        L{m} = D4./(E2-ones(1,1,Length)).*[-ones(1,1,Length)-E2 2*E;-2*E ones(1,1,Length)+E2]; % mth layer stiffness matrix [3], Eq. (23)
    end
end
M{1} = L{1};
if  strcmp(MatrixMethod,'TMM') % layer transfer matrix multiplication routine, starting with the uppermost layer to obtain the laminate transfer matrix [5], Eq. (50)
    for m = 2:SuperLayerSize
        M{1} = pagemtimes(M{1},L{m});
    end
    for m = 1:length(Pattern)
        M{m+1} = pagemtimes(M{m},M{Pattern(m)});
    end
    if  SymmetricSystem
        M2{1} = L{end};
        for m = SuperLayerSize-1:-1:1
            M2{1} = pagemtimes(M2{1},L{m});
        end
        for m = 1:length(Pattern)
            M2{m+1} = pagemtimes(M2{m},M2{Pattern(m)});
        end
        M{end} = pagemtimes(M{end},M2{end});
    end
    Y = reshape(M{end}(2,1,:),Size); % dispersion equation [2], Eq. (11)
elseif strcmp(MatrixMethod,'SMM') % recursive algorithm to obtain the laminate stiffness matrix [4], here in the 2x2 form for shear horizontal motion
    for m = 2:SuperLayerSize % computing the super layer stiffness matrix using the recursive algorithm [4], Eq. (21)
        M0 = L{m}(1,1,:)-M{1}(2,2,:);
        M1 = M{1}(1,2,:)./M0;
        M2 = L{m}(2,1,:)./M0;
        M{1} = [M{1}(1,1,:)+M1.*M{1}(2,1,:) -M1.*L{m}(1,2,:);M2.*M{1}(2,1,:) L{m}(2,2,:)-M2.*L{m}(1,2,:)];
    end
    for m = 1:length(Pattern) % repeating [4], Eq. (21) on the super layer stiffness matrix in an efficient manner (less loop iterations than super layer repetitions contained in the laminate)
        M0 = M{Pattern(m)}(1,1,:)-M{m}(2,2,:);
        M1 = M{m}(1,2,:)./M0;
        M2 = M{Pattern(m)}(2,1,:)./M0;
        M{m+1} = [M{m}(1,1,:)+M1.*M{m}(2,1,:) -M1.*M{Pattern(m)}(1,2,:);M2.*M{m}(2,1,:) M{Pattern(m)}(2,2,:)-M2.*M{Pattern(m)}(1,2,:)];
    end
    if  SymmetricSystem % calculating the laminate stiffness matrix from the upper half of the layup [5], Eq. (24); this is an alternative to [4], Eq. (43)
        M1 = -M{end}(1,2,:)./(2*M{end}(2,2,:));
        M{end} = [M{end}(1,1,:)+M1.*M{end}(2,1,:) M1.*M{end}(2,1,:);-M1.*M{end}(2,1,:) -M{end}(1,1,:)-M1.*M{end}(2,1,:)];
    end
    for j = 1:Length
        Y(j) = det(M{end}(:,:,j)); % dispersion equation [4], Eq. (38)
    end
end
Y = abs(Y); % absolute value of the dispersion equation amplitude
toc

%% plot dispersion diagram
f = figure('Name','Dispersion diagram','Units','normalized','Toolbar','none','OuterPosition',[0 0 1 1],'color','w');
imagesc(Frequency,PhaseVelocity/1e3,20*log10(Y)) % dispersion equation amplitude in decibels
colormap(gray(1024));
ax = gca;
ax.FontSize = 24;
ax.Title.FontSize = 30;
ax.XLabel.FontSize = 30;
ax.YLabel.FontSize = 30;
ax.Title.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
if  SymmetricSystem
    if  Repetitions == 1
        ax.Title.String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',replace(Name,'_','\_'),' [',replace(num2str(LayerOrientations),whitespacePattern,'/'),']$_{\mathrm s}$ using ',MatrixMethod];
    else
        ax.Title.String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',replace(Name,'_','\_'),' [',replace(num2str(LayerOrientations),whitespacePattern,'/'),']$_{',num2str(Repetitions),'\mathrm s}$ using ',MatrixMethod];
    end
else
    if  Repetitions == 1
        ax.Title.String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',replace(Name,'_','\_'),' [',replace(num2str(LayerOrientations),whitespacePattern,'/'),'] using ',MatrixMethod];
    else
        ax.Title.String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',replace(Name,'_','\_'),' [',replace(num2str(LayerOrientations),whitespacePattern,'/'),']$_{',num2str(Repetitions),'}$ using ',MatrixMethod];
    end
end
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.YDir = 'normal';
if  strcmp(MatrixMethod,'TMM')
    ax.CLim(2) = 2.5*ax.CLim(1); % compensate the exponentially growing dispersion equation amplitudes caused by the numerical instability
elseif strcmp(MatrixMethod,'SMM')
    ax.CLim = ax.CLim.*[1.15 .9];
end
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
datacursormode on
d.Interpreter = 'latex';
d.UpdateFcn = @MyCursor;

function output_txt = MyCursor(~,event_obj)
    output_txt = {['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
end