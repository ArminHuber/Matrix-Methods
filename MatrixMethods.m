%% Matrix Methods
% Update 2026/05/23
% written by: Armin Huber, armin.huber@dlr.de
% -------------------------------------------------------------------------
% Calculate guided wave dispersion diagrams for free anisotropic
% multilayered plates by using the transfer matrix method (TMM) or the
% stiffness matrix method (SMM). The dispersion equation amplitude is
% scaled to a grayscale. Modal solutions are represented by minima, i.e.,
% by dark shading. Notice in the dispersion diagram how TMM becomes
% numerically unstable with increasing frequency.

% References
% [1] A. H. Nayfeh, "The general problem of elastic wave propagation in
%     multilayered anisotropic media," J. Acoust. Soc. Am. 89(4), 1521–1531
%     (1991).
% [2] S. I. Rokhlin and L. Wang, "Stable recursive algorithm for elastic
%     wave propagation in layered anisotropic media: Stiffness matrix
%     method," J. Acoust. Soc. Am. 112(3), 822-834 (2002).
% [3] A. M. A. Huber, "Classification of solutions for guided waves in
%     fluid-loaded viscoelastic composites with large numbers of layers,"
%     J. Acoust. Soc. Am. 154(2), 1073–1094 (2023).
% [4] Guo et al., "Guided waves propagation in arbitrarily stacked 
%     composite laminates: Between-layers incompatibility issue resolution 
%     using hybrid matrix strategy," Compos. Struct. 322 (2023).

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

LayerOrientations = [0 90 45 -45]; % orientation of the layers in the super layer (deg), e.g., [0 90 -45 45]
LayerThicknesses = [.125 .125 .125 .125]; % corresponding layer thicknesses (mm), e.g., [.125 .125 .125 .125]

Repetitions = 2; % number of repetitions of the super layer, e.g., 2 for [0 90]2 = [0 90 0 90]
SymmetricSystem = 1; % 1: yes 0: no; if yes, the laminate is mirrored, e.g., [0 90]2s = [0 90 0 90 90 0 90 0]

MatrixMethod = 'TMM'; % select matrix method: 'TMM' or 'SMM' (transfer matrix method or stiffness matrix method)

PropagationAngle = 0; % wave propagation direction with respect to LayerOrientations (deg)
PhaseVelocityLimit = 20; % (m/ms)
FrequencyLimit = 2000; % (kHz)
PhaseVelocitySteps = 500;
FrequencySteps = 500;

%% calculate transformed stiffnesses
C = C*1e9; % stiffness (Pa)
SuperLayerSize = length(LayerOrientations);
Beta = LayerOrientations-PropagationAngle; % angle between the fibers of each layer and the propagation angle
for m = 1:SuperLayerSize % transformed stiffnesses for orthotropic materials for rotation about x3 [1], Appendix A
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
    if  c{m}(1,6) == 0 % for propagation along or normal to the fibers, c16, c26, c36, c45 remain zero; to avoid infinities and NaNs in the below computation, we give them a small value; it is not elegant but it works; a clean solution is given in [4]
        c{m}(1,6) = 1;
        c{m}(2,6) = 1;
        c{m}(3,6) = 1;
        c{m}(4,5) = 1;
    end
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
Y = complex(zeros(Size)); % preallocate memory for the dispersion equation amplitude
k = reshape(Wavenumber,1,1,Length); % reshape into pages in order to use the pagewise matrix functions "pagemtimes" and "pagemrdivide"; these enable implicit parallel computing and we avoid slow loops
k2 = k.^2;
k4 = k2.^2;
k6 = k2.^3;
rw2 = reshape(repmat(Density*AngularFrequency.^2,Size(1),1),1,1,Length);
r2w4 = rw2.^2;
r3w6 = rw2.^3;
for m = 1:SuperLayerSize
    Delta = c{m}(3,3)*c{m}(4,4)*c{m}(5,5)-c{m}(3,3)*c{m}(4,5)^2; % the frequency and phase velocity independent components of the polynomial coefficients A1, A2, A3 in [3], Eq. (11) and Appendix B
    a11 = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta;
    a12 = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta;
    a21 = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta;
    a22 = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta;
    a23 = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta;
    a31 = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta;
    a32 = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta;
    a33 = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta;
    a34 = -1/Delta;
    A1 = a11*k2+a12*rw2; % polynomial coefficients in [3], Eq. (11) and Appendix B
    A2 = a21*k4+a22*rw2.*k2+a23*r2w4;
    A3 = a31*k6+a32*rw2.*k4+a33*r2w4.*k2+a34*r3w6;
    d1 = A1/3;
    d2 = A2/3-d1.^2;
    d3 = d1.^3-d1.*A2/2+A3/2;
    d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
    d5 = d2./d4;
    d6 = (d5-d4)/2-d1;
    d7 = (d5+d4)/2i*sqrt(3);
    k32 = [d6+d7 d6-d7 d4-d5-d1];
    k3 = sqrt(k32); % out-of-plane wavenumber components of the three pairs of bulk waves (quasi-L+-, quasi-SV+-, quasi-SH+-) in [3], Eq. (12)
    k3k = k3.*k;
    m11 = c{m}(1,1)*k2+c{m}(5,5)*k32-rw2; % components of the Christoffel matrix [3], Eq. (10)
    m22 = c{m}(6,6)*k2+c{m}(4,4)*k32-rw2;
    m12 = c{m}(1,6)*k2+c{m}(4,5)*k32;
    m13 = (c{m}(1,3)+c{m}(5,5))*k3k;
    m23 = (c{m}(3,6)+c{m}(4,5))*k3k;
    m1 = m13.*m22-m12.*m23;
    V = (m11.*m23-m13.*m12)./m1; % displacement amplitude (polarization) in the x2 direction [3], Eq. (13)
    W = (m11.*m22-m12.^2)./-m1; % displacement amplitude (polarization) in the x3 direction [3], Eq. (14)
    e1 = k.*W+k3;
    e2 = k3.*V;
    D3 = 1i*((c{m}(1,3)+c{m}(3,6)*V).*k+c{m}(3,3)*k3.*W); % out-of-plane stress amplitude (sigma33) [3], Eq. (15)
    D4 = 1i*(c{m}(4,5)*e1+c{m}(4,4)*e2); % shear stress amplitude (sigma23) [3], Eq. (15)
    D5 = 1i*(c{m}(5,5)*e1+c{m}(4,5)*e2); % shear stress amplitude (sigma13) [3], Eq. (15)
    Phi = 1i*k3*LayerThicknesses(m);
    E = exp(Phi); % harmonic term in the out-of-plane direction in [3], Eqs. (16)-(17)
    if  strcmp(MatrixMethod,'TMM')
        E_ = exp(-Phi);
        L1 = [E E_;V.*E V.*E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_;D4.*E -D4.*E_]; % first 6x6 matrix in [3], Eq. (20)
        L2 = [ones(1,6,Length);V V;W -W;D3 D3;D5 -D5;D4 -D4]; % second 6x6 matrix in [3], Eq. (20)
    elseif strcmp(MatrixMethod,'SMM')
        L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4]; % first 6x6 matrix in [3], Eq. (21)
        L2 = [ones(1,3,Length) E;V V.*E;W -W.*E;E ones(1,3,Length);V.*E V;W.*E -W]; % second 6x6 matrix in [3], Eq. (21)
    end
    L{m} = pagemrdivide(L1,L2); % mth layer transfer/stiffness matrix [3], Eqs. (20)/(21)
end
M{1} = L{1};
if  strcmp(MatrixMethod,'TMM') % layer transfer matrix multiplication routine, starting with the uppermost layer to obtain the laminate transfer matrix [3], Eq. (50)
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
    for j = 1:Length
        Y(j) = det(M{end}(4:6,1:3,j)); % dispersion equation [1], Eq. (27); there is no "pagedet" function yet, so we need to loop
    end
elseif strcmp(MatrixMethod,'SMM') % recursive algorithm to obtain the laminate stiffness matrix [2]
    for m = 2:SuperLayerSize % computing the super layer stiffness matrix using the recursive algorithm [2], Eq. (21)
        M0 = L{m}(1:3,1:3,:)-M{1}(4:6,4:6,:);
        M1 = pagemrdivide(M{1}(1:3,4:6,:),M0);
        M2 = pagemrdivide(L{m}(4:6,1:3,:),M0);
        M{1} = [M{1}(1:3,1:3,:)+pagemtimes(M1,M{1}(4:6,1:3,:)) -pagemtimes(M1,L{m}(1:3,4:6,:));pagemtimes(M2,M{1}(4:6,1:3,:)) L{m}(4:6,4:6,:)-pagemtimes(M2,L{m}(1:3,4:6,:))];
    end
    for m = 1:length(Pattern) % repeating [2], Eq. (21) on the super layer stiffness matrix in an efficient manner (less loop iterations than super layer repetitions contained in the laminate)
        M0 = M{Pattern(m)}(1:3,1:3,:)-M{m}(4:6,4:6,:);
        M1 = pagemrdivide(M{m}(1:3,4:6,:),M0);
        M2 = pagemrdivide(M{Pattern(m)}(4:6,1:3,:),M0);
        M{m+1} = [M{m}(1:3,1:3,:)+pagemtimes(M1,M{m}(4:6,1:3,:)) -pagemtimes(M1,M{Pattern(m)}(1:3,4:6,:));pagemtimes(M2,M{m}(4:6,1:3,:)) M{Pattern(m)}(4:6,4:6,:)-pagemtimes(M2,M{Pattern(m)}(1:3,4:6,:))];
    end
    if  SymmetricSystem % calculating the laminate stiffness matrix from the upper half of the layup [3], Eq. (24); this is an alternative to [2], Eq. (43)
        I = [1 1 -1;-1 -1 1;-1 -1 1];
        M0 = M{end}(4:6,4:6,:).*I-M{end}(4:6,4:6,:);
        M1 = pagemrdivide(M{end}(1:3,4:6,:),M0);
        M2 = pagemrdivide(M{end}(1:3,4:6,:).*I,M0);
        M{end} = [M{end}(1:3,1:3,:)+pagemtimes(M1,M{end}(4:6,1:3,:)) -pagemtimes(M1,M{end}(4:6,1:3,:).*I);pagemtimes(M2,M{end}(4:6,1:3,:)) M{end}(1:3,1:3,:).*I-pagemtimes(M2,M{end}(4:6,1:3,:).*I)];
    end
    for j = 1:Length
        Y(j) = det(M{end}(:,:,j)); % dispersion equation [2], Eq. (38)
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