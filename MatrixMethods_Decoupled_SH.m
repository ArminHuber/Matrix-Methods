%% Matrix Methods Decoupled SH
% Update 2023/12/13
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

LayerOrientations = [0 90]; % orientation of the layers in the unit cell (deg), e.g., [0 90]
LayerThicknesses = [.125 .125]; % corresponding layer thicknesses (mm), e.g., [.125 .125]

Repetitions = 2; % number of repetitions of the unit cell, e.g., 2 for [0 90]2 = [0 90 0 90]
Symmetric = 1; % 1: yes 0: no; if yes, the laminate is mirrored, e.g., [0 90]2s = [0 90 0 90 90 0 90 0]

Method = 'TMM'; % select matrix method: 'TMM' or 'SMM' (transfer matrix method or stiffness matrix method)

PropagationAngle = 0; % wave propagation direction with respect to LayerOrientations (deg)
PhaseVelocityLimit = 20; % (m/ms)
FrequencyLimit = 4000; % (kHz)
PhaseVelocitySteps = 500;
FrequencySteps = 500;

%% calculate transformed stiffnesses
C = C*1e9; % stiffness (Pa)
UnitCellSize = length(LayerOrientations); %#ok<*NBRAK2> 
Beta = LayerOrientations-PropagationAngle; % angle between the fibers of each layer and the propagation angle
if  any(mod(Beta,90)) % check if the layup supports decoupling and wave propagation is along an axis of symmetry
    errordlg('Unsupported layup or propagation direction. "LayerOrientations" may only contain single layers or cross plies such as 0 and 90 degree layers and "PropagationAngle" may only be along or normal to a principal axis, e.g., 0 or 90 degrees in this case.','Error');
    return
end
for m = 1:UnitCellSize % transformed stiffnesses for orthotropic materials for rotation about x3 [1], Appendix A
    s = sind(Beta(m));
    g = cosd(Beta(m));
    c{m}(4,4) = C(4,4)*g^2+C(5,5)*s^2; %#ok<*SAGROW> 
    c{m}(6,6) = C(6,6)+(C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2;
end
LayerThicknesses = LayerThicknesses/1e3; % (m)
Layers = Repetitions*UnitCellSize;
Thickness = Repetitions*sum(LayerThicknesses)*1e3;
if  Symmetric %#ok<*UNRCH>
    Layers = 2*Layers; 
    Thickness = 2*Thickness;
end
disp(['Method:    ',Method,newline...
      'Layers:    ',num2str(Layers),newline...
      'Thickness: ',num2str(Thickness),' mm'])
warning('off','MATLAB:nearlySingularMatrix') % turn warnings off
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:illConditionedMatrix')

%% computation
tic
Frequency = 0:FrequencyLimit/FrequencySteps:FrequencyLimit; Frequency(1) = 1e-3; % generate a range of frequencies
PhaseVelocity = (0:PhaseVelocityLimit/PhaseVelocitySteps:PhaseVelocityLimit)*1e3; PhaseVelocity(1) = 1e-3; % generate a range of phase velocities
AngularFrequency = 2*pi*Frequency*1e3;
Wavenumber = AngularFrequency./PhaseVelocity';
I = -1; % used for symmetric layups in SMM to obtain the global stiffness matrix from the matrix of the half layup
Y = zeros(length(PhaseVelocity),length(Frequency)); % preallocate memory
h = waitbar(0,sprintf('0 of %d (0 %%)',length(Frequency)),'Name','Calculating frequency...');
for i = 1:length(Frequency)
    for j = 1:length(PhaseVelocity)
        for m = 1:UnitCellSize
            Alpha = sqrt((Density*PhaseVelocity(j)^2-c{m}(6,6))/c{m}(4,4)); % alpha from [2], Eq. (5); alpha gives the out-of-plane wavenumber component ratios for the pair of bulk waves (SH+-), i.e., their propagation directions in the sagittal plane x1-x3
            D4 = Alpha*c{m}(4,4); % shear stress amplitude (sigma23) [2], Eq. (8)
            if  strcmp(Method,'TMM')
                G = Wavenumber(j,i)*Alpha*LayerThicknesses(m);
                L{m} = [cos(G) 1i*sin(G)/D4;1i*D4*sin(G) cos(G)]; % mth layer transfer matrix [2], Eq. (7), relating the displacements and stresses at the top of the mth layer to those at its bottom
            elseif strcmp(Method,'SMM')
                H = exp(1i*Wavenumber(j,i)*Alpha*LayerThicknesses(m)); % harmonic term in the out-of-plane direction [1], Eq. (15)
                L{m} = 1i*Wavenumber(j,i)*D4/(H^2-1)*[-1-H^2 2*H;-2*H 1+H^2]; % mth layer stiffness matrix [3], Eq. (23), relating the stresses at the top and bottom of the mth layer to the corresponding displacements
            end
        end
        if  strcmp(Method,'TMM') % layer transfer matrix multiplication routine, starting with the lowermost layer to obtain the laminate transfer matrix [1], Eq. (18)
            if  Symmetric
                M = L{1};
                for m = 2:UnitCellSize
                    M = M*L{m};
                end
                for n = 2:Repetitions
                    for m = 1:UnitCellSize
                        M = M*L{m};
                    end
                end
                for n = 1:Repetitions
                    for m = UnitCellSize:-1:1
                        M = M*L{m};
                    end
                end
            else
                if  UnitCellSize == 1
                    M = L{1}^Repetitions;
                else
                    M = L{end};
                    for m = UnitCellSize-1:-1:1
                        M = M*L{m};
                    end
                    for n = 2:Repetitions
                        for m = UnitCellSize:-1:1
                            M = M*L{m};
                        end
                    end
                end
            end
            Y(j,i) = abs(det(M(2,1))); % dispersion equation (absolute value) [2], Eq. (11)
        elseif strcmp(Method,'SMM') % recursive algorithm to obtain the laminate stiffnes matrix [4], here in the 2x2 form for shear horizontal motion
            M = L{1};
            for m = 2:UnitCellSize % computation of the unit cell stiffness matrix using the recursive algorithm [4], Eq. (21)
                N = inv(L{m}(1,1)-M(2,2));
                M = [M(1,1)+M(1,2)*N*M(2,1) -M(1,2)*N*L{m}(1,2);L{m}(2,1)*N*M(2,1) L{m}(2,2)-L{m}(2,1)*N*L{m}(1,2)]; %#ok<*MINV> 
            end
            MM{1} = M;
            for m = 2:Repetitions % repeating [4], Eq. (21) on the unit cell stiffness matrix as many times as repetitions are contained in the laminate
                N = inv(MM{1}(1,1)-MM{m-1}(2,2));
                MM{m} = [MM{m-1}(1,1)+MM{m-1}(1,2)*N*MM{m-1}(2,1) -MM{m-1}(1,2)*N*MM{1}(1,2);MM{1}(2,1)*N*MM{m-1}(2,1) MM{1}(2,2)-MM{1}(2,1)*N*MM{1}(1,2)];
            end
            if  Symmetric % calculating the laminate stiffness matrix from the upper half of the layup [5], Eq. (24); this is an alternative to [4], Eq. (43)
                N = inv(MM{end}(2,2)*I-MM{end}(2,2));
                MM{end} = [MM{end}(1,1)+MM{end}(1,2)*N*MM{end}(2,1) -MM{end}(1,2)*N*(MM{end}(2,1)*I);(MM{end}(1,2)*I)*N*MM{end}(2,1) (MM{end}(1,1)*I)-(MM{end}(1,2)*I)*N*(MM{end}(2,1)*I)];
            end
            Y(j,i) = abs(det(MM{end})); % dispersion equation (absolute value) [4], Eq. (38)
        end
    end
    waitbar(i/length(Frequency),h,sprintf('%d of %d (%.0f %%), elapsed %.0f sec',i,length(Frequency),100*i/length(Frequency),toc))
end
close(h)
toc

%% dispersion diagram
f = figure('Name','Dispersion diagram','Units','normalized','Toolbar','none','OuterPosition',[0 0 1 1],'color','w');
datacursormode on
imagesc(Frequency,PhaseVelocity/1e3,20*log10(Y)) % dispersion equation amplitude in decibel
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
if  ~Symmetric
    if  Repetitions == 1
        ax.Title.String = ['SH mode dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',char(join(split(Name,'_'),'\_')),' [',char(join(split(num2str(LayerOrientations)),'/')),'] using ',Method];
    else
        ax.Title.String = ['SH mode dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',char(join(split(Name,'_'),'\_')),' [',char(join(split(num2str(LayerOrientations)),'/')),']$_{',num2str(Repetitions),'}$ using ',Method];
    end
else
    if  Repetitions == 1
        ax.Title.String = ['SH mode dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',char(join(split(Name,'_'),'\_')),' [',char(join(split(num2str(LayerOrientations)),'/')),']$_{\mathrm s}$ using ',Method];
    else
        ax.Title.String = ['SH mode dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',char(join(split(Name,'_'),'\_')),' [',char(join(split(num2str(LayerOrientations)),'/')),']$_{',num2str(Repetitions),'\mathrm s}$ using ',Method];
    end
end
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.YDir = 'normal';
if  strcmp(Method,'TMM')
    ax.CLim(2) = 8*ax.CLim(1); % compensate the exponentially growing dispersion equation amplitudes caused by the numerical instability
end
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @MyCursor;

function output_txt = MyCursor(~,event_obj)
    output_txt = {['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
end