%% Matrix Methods
% Update 2023/12/13
% written by: Armin Huber, armin.huber@dlr.de
% -------------------------------------------------------------------------
% Calculate guided wave dispersion diagrams for free anisotropic
% multilayered plates by using the transfer matrix method (TMM) or the
% stiffness matrix method (SMM). The dispersion equation amplitude is
% scaled to a grayscale. Modal solutions are represented by minima, i.e.,
% by dark shading. Notice in the dispersion diagram how TMM becomes
% numerically unstable with increasing frequency.

% [1] A. H. Nayfeh, "The general problem of elastic wave propagation in
%     multilayered anisotropic media," J. Acoust. Soc. Am. 89(4), 1521–1531
%     (1991).
% [2] S. I. Rokhlin and L. Wang, "Stable recursive algorithm for elastic
%     wave propagation in layered anisotropic media: Stiffness matrix
%     method," J. Acoust. Soc. Am. 112(3), 822-834 (2002).
% [3] A. M. A. Huber, "Classification of solutions for guided waves in
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

LayerOrientations = [0 90]; % orientation of the layers in the unit cell (deg), e.g., [0 90 -45 45]
LayerThicknesses = [.125 .125]; % corresponding layer thicknesses (mm), e.g., [.125 .125 .125 .125]

Repetitions = 2; % number of repetitions of the unit cell, e.g., 2 for [0 90]2 = [0 90 0 90]
Symmetric = 1; % 1: yes 0: no; if yes, the laminate is mirrored, e.g., [0 90]2s = [0 90 0 90 90 0 90 0]

Method = 'TMM'; % select matrix method: 'TMM' or 'SMM' (transfer matrix method or stiffness matrix method)

PropagationAngle = 30; % wave propagation direction with respect to LayerOrientations (deg)
PhaseVelocityLimit = 20; % (m/ms)
FrequencyLimit = 4000; % (kHz)
PhaseVelocitySteps = 500;
FrequencySteps = 500;

%% calculate transformed stiffnesses
C = C*1e9; % stiffness (Pa)
UnitCellSize = length(LayerOrientations); %#ok<*NBRAK2> 
Beta = LayerOrientations-PropagationAngle; % angle between the fibers of each layer and the propagation angle
for m = 1:UnitCellSize % transformed stiffnesses for orthotropic materials for rotation about x3 [1], Appendix A
    s = sind(Beta(m));
    g = cosd(Beta(m));
    c{m}(1,1) = C(1,1)*g^4+C(2,2)*s^4+2*(C(1,2)+2*C(6,6))*s^2*g^2; %#ok<*SAGROW> 
    c{m}(1,2) = (C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2+C(1,2);
    c{m}(1,3) = C(1,3)*g^2+C(2,3)*s^2;
    c{m}(1,6) = (C(1,2)+2*C(6,6)-C(1,1))*s*g^3+(C(2,2)-C(1,2)-2*C(6,6))*g*s^3;
    c{m}(2,2) = C(1,1)*s^4+C(2,2)*g^4+2*(C(1,2)+2*C(6,6))*s^2*g^2;
    c{m}(2,3) = C(2,3)*g^2+C(1,3)*s^2;
    c{m}(2,6) = (C(1,2)+2*C(6,6)-C(1,1))*g*s^3+(C(2,2)-C(1,2)-2*C(6,6))*s*g^3;
    c{m}(3,3) = C(3,3);
    c{m}(3,6) = (C(2,3)-C(1,3))*s*g;
    c{m}(4,4) = C(4,4)*g^2+C(5,5)*s^2;
    c{m}(4,5) = (C(4,4)-C(5,5))*s*g;
    c{m}(5,5) = C(5,5)*g^2+C(4,4)*s^2;
    c{m}(6,6) = C(6,6)+(C(1,1)+C(2,2)-2*C(1,2)-4*C(6,6))*s^2*g^2;
    if  c{m}(1,6) == 0 % for propagation along or normal to the fibers, c16, c26, c36, c45 remain zero; to avoid infinities and NANs in the below computation, we give them a small value; it is not elegant but it works
        c{m}(1,6) = 1;
        c{m}(2,6) = 1;
        c{m}(3,6) = 1;
        c{m}(4,5) = 1;
    end
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
for m = 1:UnitCellSize % calculating the frequency and phase velocity independent components of the polynomial coefficients A1, A2, A3 in [1], Eq. (10) and Appendix B; notice that there is one incorrect sign in A2 in [1], Appendix B; it is corrected in a21 below
    Delta(m) = c{m}(3,3)*c{m}(4,4)*c{m}(5,5)-c{m}(3,3)*c{m}(4,5)^2;
    a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
    a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
    a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
    a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
    a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
    a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
    a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
    a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
    a34(m) = -1/Delta(m);
end
I = [1 1 -1;-1 -1 1;-1 -1 1]; % used for symmetric layups in SMM to obtain the global stiffness matrix from the matrix of the half layup [3], Eq. (25)
Y = zeros(length(PhaseVelocity),length(Frequency)); % preallocate memory
h = waitbar(0,sprintf('0 of %d (0 %%)',length(Frequency)),'Name','Calculating frequency...');
for i = 1:length(Frequency)
    for j = 1:length(PhaseVelocity)
        for m = 1:UnitCellSize
            rc2 = Density*PhaseVelocity(j)^2;
            r2c4 = rc2^2;
            A1 = a11(m)+a12(m)*rc2; % polynomial coefficients in [1], Eq. (10) and Appendix B
            A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
            A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*rc2^3;
            b1 = A2/3-A1^2/9;
            b2 = A1^3/27-A1*A2/6+A3/2;
            b3 = (sqrt(b2^2+b1^3)-b2)^(1/3);
            b4 = b1/(2*b3)-b3/2;
            b5 = b1/b3;
            b6 = (sqrt(3)*(b3+b5)*1i)/2;
            Alpha(1) = sqrt(b4-b6-A1/3); % alpha from [1], Eq. (11); alpha gives the out-of-plane wavenumber component ratios for the three pairs of bulk waves (quasi-L+-, quasi-SV+-, quasi-SH+-), i.e., their propagation directions in the sagittal plane x1-x3
            Alpha(2) = sqrt(b4+b6-A1/3);
            Alpha(3) = -sqrt(b3-b5-A1/3);
            Alpha2 = Alpha.^2;
            m11 = c{m}(1,1)-rc2+c{m}(5,5)*Alpha2; % components of the Christoffel matrix [1], Eq. (9b)
            m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
            m13 = (c{m}(1,3)+c{m}(5,5))*Alpha;
            m22 = c{m}(6,6)-rc2+c{m}(4,4)*Alpha2;
            m23 = (c{m}(3,6)+c{m}(4,5))*Alpha;
            V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23); % displacement amplitude (polarization) in the x2 direction [1], Eq. (12)
            W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13); % displacement amplitude (polarization) in the x3 direction [1], Eq. (13); W is derived differently here, but is equivalent to [1], Eq. (13)
            D3 = c{m}(1,3)+c{m}(3,6)*V+c{m}(3,3)*Alpha.*W; % stress amplitude (sigma33) in x3 (out-of-plane) [1], Eq. (15); we omit the common factor i*xi here
            D4 = c{m}(4,5)*(Alpha+W)+c{m}(4,4)*Alpha.*V; % shear stress amplitude (sigma23)
            D5 = c{m}(5,5)*(Alpha+W)+c{m}(4,5)*Alpha.*V; % shear stress amplitude (sigma13)
            if  strcmp(Method,'TMM')
                D = diag([exp(1i*Wavenumber(j,i)*[Alpha -Alpha]*LayerThicknesses(m))]); % harmonic term diagonal matrix
                X = [1 1 1 1 1 1;V V;W -W;D3 D3;D5 -D5;D4 -D4]; % displacement and stress matrix from [1], Eq. (14), but with different order in the columns for convenience
                L{m} = X*D/X; % mth layer transfer matrix [1], Eq. (17), relating the displacements and stresses at the top of the mth layer to those at its bottom
            elseif strcmp(Method,'SMM')
                H = exp(1i*Wavenumber(j,i)*Alpha*LayerThicknesses(m)); % harmonic term in the out-of-plane direction [1], Eq. (15)
                P = [1 1 1 H;V V.*H;W -W.*H;H 1 1 1;V.*H V;W.*H -W]; % displacement amplitude matrix [2], Eq. (3) and [3], Appendix C
                D = [D3 D3.*H;D5 -D5.*H;D4 -D4.*H;D3.*H D3;D5.*H -D5;D4.*H -D4]; % stress amplitude matrix [2], Eq. (5) and [3], Appendix C
                L{m} = D/P; % mth layer stiffness matrix [2], Eqs. (6)-(7), relating the stresses at the top and bottom of the mth layer to the corresponding displacements
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
            Y(j,i) = abs(det(M(4:6,1:3))); % dispersion equation (absolute value) [1], Eq. (27)
        elseif strcmp(Method,'SMM') % recursive algorithm to obtain the laminate stiffnes matrix [2]
            M = L{1};
            for m = 2:UnitCellSize % computation of the unit cell stiffness matrix using the recursive algorithm [2], Eq. (21)
                N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)]; %#ok<*MINV> 
            end
            MM{1} = M;
            for m = 2:Repetitions % repeating [2], Eq. (21) on the unit cell stiffness matrix as many times as repetitions are contained in the laminate
                N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
            end
            if  Symmetric % calculating the laminate stiffness matrix from the upper half of the layup [3], Eq. (24); this is an alternative to [2], Eq. (43)
                N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
            end
            Y(j,i) = abs(det(MM{end})); % dispersion equation (absolute value) [2], Eq. (38)
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
        ax.Title.String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',char(join(split(Name,'_'),'\_')),' [',char(join(split(num2str(LayerOrientations)),'/')),'] using ',Method];
    else
        ax.Title.String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',char(join(split(Name,'_'),'\_')),' [',char(join(split(num2str(LayerOrientations)),'/')),']$_{',num2str(Repetitions),'}$ using ',Method];
    end
else
    if  Repetitions == 1
        ax.Title.String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',char(join(split(Name,'_'),'\_')),' [',char(join(split(num2str(LayerOrientations)),'/')),']$_{\mathrm s}$ using ',Method];
    else
        ax.Title.String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness),'\,mm ',char(join(split(Name,'_'),'\_')),' [',char(join(split(num2str(LayerOrientations)),'/')),']$_{',num2str(Repetitions),'\mathrm s}$ using ',Method];
    end
end
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.YDir = 'normal';
if  strcmp(Method,'TMM')
    ax.CLim(2) = 4*ax.CLim(1); % compensate the exponentially growing dispersion equation amplitudes caused by the numerical instability
end
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @MyCursor;

function output_txt = MyCursor(~,event_obj)
    output_txt = {['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{p}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
end