% Script: epsif_calc.m
% 
% This file is supplemental to Hantsoo et al., 2023 (JGR: Biogeosciences).
% script written by K. Hantsoo, Nov. 2022.
%
% This script computes S8-pyrite isotopic offset using given input values q
% and n. It assumes that pyrite forms by the polysulfide reaction detailed
% by Luther (1991) and Butler et al. (2004).
% 
% PSIF calculations for polysulfide chain lengths S4 through S7 use d34S
% values from Amrani et al. (2006), who measured the d34S of polsyulfide
% chains sorted by length and obtained distinct d34S values for S4, S5, S6,
% and S7. PSIF calcs for S9 and its disproportionation to S8 use the 3.4 to
% 5.4 per mil fractionation observed by Amrani and Aizenshtat (2004).


% References:
%
% Amrani, A., & Aizenshtat, Z. (2004). Mechanisms of sulfur introduction
% chemically controlled: ?34S imprint. Organic Geochemistry, 35(11-12),
% 1319-1336. https://doi.org/10.1016/j.orggeochem.2004.06.019
%
% Amrani, A., Kamyshny, A., Lev, O., & Aizenshtat, Z. (2006). Sulfur stable isotope
% distribution of polysulfide anions in an (NH4)2Sn aqueous solution.
% Inorganic Chemistry, 45(4), 1427-1429. https://doi.org/10.1021/ic051748r
% 
% Butler, I. B., Böttcher, M. E., Rickard, D., & Oldroyd, A. (2004). Sulfur
% isotope partitioning during experimental formation of pyrite via the
% polysulfide and hydrogen sulfide pathways: implications for the
% interpretation of sedimentary and hydrothermal pyrite isotope records.
% Earth and Planetary Science Letters, 228(3-4), 495-509.
% https://doi.org/10.1016/j.epsl.2004.10.005
% 
% Luther III, G. W. (1991). Pyrite synthesis via polysulfide compounds.
% Geochimica et Cosmochimica Acta, 55(10), 2839-2849.
% https://doi.org/10.1016/0016-7037(91)90449-f


% Position notation is as follows:
% 
% S42-:
% S - S - S - S
% a   b   b   a
% 
% S52-:
% S - S - S - S - S
% a   b   c   b   a
% 
% S62-:
% S - S - S - S - S - S
% a   b   c   c   b   a
% 
% S72-:
% S - S - S - S - S - S - S
% a   b   c   d   c   b   a
% 
% S82-:
% S - S - S - S - S - S - S - S
% a   b   c   d   d   c   b   a
% 
% S92-:
% S - S - S - S - S - S - S - S - S
% a   b   c   d   e   d   c   b   a
% 
% 
% System of equations for computing d34S of individual d34S values by position:
% 
% d34S_S4-a = q
% d34S_S4-b = q + eb
% 
% d34S_S5-a = q + p
% d34S_S5-b = q + p + eb
% d34S_S5-c = q + p + fc*eb
% 
% d34S_S6-a = q + 2p
% d34S_S6-b = q + 2p + eb
% d34S_S6-c = q + 2p + fc*eb
% 
% d34S_S7-a = q + 3p
% d34S_S7-b = q + 3p + eb
% d34S_S7-c = q + 3p + fc*eb
% d34S_S7-d = q + 3p + fd*eb
% 
% d34S_S8-a = q + 4p
% d34S_S8-b = q + 4p + eb
% d34S_S8-c = q + 4p + fc*eb
% d34S_S8-d = q + 4p + fd*eb
% 
% d34S_S9-a = q + 5p
% d34S_S9-b = q + 5p + eb
% d34S_S9-c = q + 5p + fc*eb
% d34S_S9-d = q + 5p + fd*eb
% d34S_S9-e = q + 5p + fe*eb
%
%
% fc, fd, and fe should not be less than 1.



clear variables
close all

% Inputs: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
system_d34S = 16.5 ; % bulk d34S of all sulfur in system from Amrani et al., 2006
extrap_d34S_HS = 13.8 ; % calculated d34S of dissolved sulfide from Amrani et al., 2006
extrap_d34S_S9 = 23.0 ; % extrapolated d34S of S92- from Amrani et al., 2006
extrap_d34S_S8 = 21.8 ; % extrapolated d34S of S82- from Amrani et al., 2006

%The reported standard error of these four measurements is +/- 0.3 per mil:
d_4_obs = 16.9 ; % 16.9 is the observed value from Amrani et al., 2006
d_5_obs = 18.1 ; % 18.1 is the observed value from Amrani et al., 2006
d_6_obs = 19.2 ; % 19.2 is the observed value from Amrani et al., 2006
d_7_obs = 20.6 ; % 20.6 is the observed value from Amrani et al., 2006

d_4_in = d_4_obs - system_d34S ;
d_5_in = d_5_obs - system_d34S ;
d_6_in = d_6_obs - system_d34S ;
d_7_in = d_7_obs - system_d34S ;

S8_HS_offset = 4.4 ; % 4.4 +/- 1.0 per mil (Amrani & Aizenshtat, 2004)

% Create matrix:
lengthq = 90;
lengthp = 100;
q = linspace((extrap_d34S_HS-system_d34S),d_4_in,lengthq) ; % input values for q (see equations at top of script)
p = linspace(0,1.2,lengthp) ; % input values for p (see equations at top of script)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Run script:

for b=1:lengthq %     varying q
    for c=1:lengthp % varying p

    % Solve equations for the coefficients eb, fc, fd, and fe:
    
        % solve equation d_4 for eb:
        eb(b) = 2*d_4_in - 2*q(b) ;

        % solve equation d_5 for fc:
        fc(b,c) = (5*d_5_in - 2*eb(b) - 5*q(b) - 5*p(c))/eb(b) ; 

        % alternate solution for fc using equation d_6:
        %fc(b,c) = (6*d_6_in - 6*q(b) - 12*p(c) - 2*eb(b))/(2*eb(b)) ;

        % solve equation d_7 for fd:
        fd(b,c) = (7*d_7_in - 7*q(b) - 21*p(c) - 2*eb(b) - ...
            2*fc(b,c)*eb(b))/eb(b) ;

        % solve S92- disproportionation equation for fe:
        fe(b,c) =  (8*S8_HS_offset - 2*eb(b) - 2*fc(b,c)*eb(b) - ...
           2*fd(b,c)*eb(b)) / eb(b) ;
    
    % Compute output polysulfide values to confirm that they match
    % observations and input values:
        
        d_4_a(b,c) = q(b) ;
        d_4_b(b,c) = q(b) + eb(b) ;
        d_4_out(b,c) = 0.25*(2*d_4_a(b,c) + 2*d_4_b(b,c)) ;
        
        
        d_5_a(b,c) = q(b) + p(c) ;
        d_5_b(b,c) = q(b) + p(c) + eb(b) ;
        d_5_c(b,c) = q(b) + p(c) + fc(b,c)*eb(b) ;
        d_5_out(b,c) = 0.2*(2*d_5_a(b,c) + 2*d_5_b(b,c) + d_5_c(b,c)) ;
        
        
        d_6_a(b,c) = q(b) + 2*p(c);
        d_6_b(b,c) = q(b) + 2*p(c) + eb(b);
        d_6_c(b,c) = q(b) + 2*p(c) + fc(b,c)*eb(b);
        d_6_out(b,c) = (1/6)*(2*d_6_a(b,c) + 2*d_6_b(b,c) + 2*d_6_c(b,c));
        
        
        d_7_a(b,c) = q(b) + 3*p(c) ;
        d_7_b(b,c) = q(b) + 3*p(c) + eb(b) ;
        d_7_c(b,c) = q(b) + 3*p(c) + fc(b,c)*eb(b) ;
        d_7_d(b,c) = q(b) + 3*p(c) + fd(b,c)*eb(b) ;
        d_7_out(b,c) = (1/7)*(2*d_7_a(b,c) + 2*d_7_b(b,c) + 2*d_7_c(b,c) + d_7_d(b,c));
        
        
        d_8_a(b,c) = q(b) + 4*p(c) ;
        d_8_b(b,c) = q(b) + 4*p(c) + eb(b) ;
        d_8_c(b,c) = q(b) + 4*p(c) + fc(b,c)*eb(b) ;
        d_8_d(b,c) = q(b) + 4*p(c) + fd(b,c)*eb(b) ;
        d_8_out(b,c) = (1/8)*(2*d_8_a(b,c) + 2*d_8_b(b,c) + 2*d_8_c(b,c) + 2*d_8_d(b,c));
        
        
        d_9_a(b,c) = q(b) + 5*p(c) ;
        d_9_b(b,c) = q(b) + 5*p(c) + eb(b) ;
        d_9_c(b,c) = q(b) + 5*p(c) + fc(b,c)*eb(b) ;
        d_9_d(b,c) = q(b) + 5*p(c) + fd(b,c)*eb(b) ;
        d_9_e(b,c) = q(b) + 5*p(c) + fe(b,c)*eb(b) ;
        d_9_out(b,c) = (1/9)*(2*d_9_a(b,c) + 2*d_9_b(b,c) + 2*d_9_c(b,c) + 2*d_9_d(b,c) + d_9_e(b,c));
        
        
        d_8_postdisprop(b,c) = 1/8 * (d_9_a(b,c) + 2*d_9_b(b,c) + ...
            2*d_9_c(b,c) + 2*d_9_d(b,c) + d_9_e(b,c));
        
        
    % Calculate d34S of pyrite derived from S52-:
        d_pyr_S5(b,c) = 0.5*(d_5_a(b,c) + d_5_b(b,c)) ;
    % Now the same for S42- and S62-:
        d_pyr_S4(b,c) = 0.5*(d_4_a(b,c) + d_4_b(b,c)) ;
        d_pyr_S6(b,c) = 0.5*(d_6_a(b,c) + d_6_b(b,c)) ;

    % per mil difference between disproportionated S8 and pyrite:
        diff_dS8_dpyr(b,c) = d_8_postdisprop(b,c) - d_pyr_S5(b,c);
        %diff_dS8_dpyr(b,c) = d_8_out(b,c) - d_pyr_S6(b,c);
    end
end

% Note: the following _diff values are the same regardless of chain length
de_diff = d_9_e - d_9_d; 
cd_diff = d_9_d - d_9_c;
bc_diff = d_9_c - d_9_b;
ab_diff = d_9_b - d_9_a;

[P,Q] = meshgrid(p,q);

figure(1)
[F,G] = contour(Q,P,diff_dS8_dpyr,[1 2 3 4 5 6 7 8 9],'k');
caxis([1 9]);
clabel(F,G)
hold on

% plot range in which computed d_9_out matches extrapolated value from Amrani dataset:
extrap_dS9_min = extrap_d34S_S9 - system_d34S - 0.3 ;
extrap_dS9_max = extrap_d34S_S9 - system_d34S + 0.3 ;
levels=[extrap_dS9_min extrap_dS9_max];
contour(Q,P,d_9_out,levels,':b')
hold on

% plot range in which eb is between 3.4 and 5.4 per mil:
ebcol = rot90(eb);
ebmatrix = repmat(ebcol,1,100);
levels1 = [3.4 5.4];
contour(Q,P,ebmatrix,levels1,'-.b')

% plot range in which d34S_e is between 0 and 5 per mil greater than d34S_d:
levels2 = [0 5];
contour(Q,P,de_diff,levels2,'-b')

% plot range in which d34S_c is between 0 and 5 per mil greater than d34S_b:
levels3 = [0 5];
contour(Q,P,bc_diff,levels3,'--b')

% plot range in which computed d34S_S6 matches measured d34S_S6
obs_dS6_min = d_6_in - 0.3 ; %2.4
obs_dS6 = d_6_in ; %=2.7
obs_dS6_max = d_6_in + 0.3 ; %3.0
levels3 = [obs_dS6_min obs_dS6 obs_dS6_max];
contour(Q,P,d_6_out,levels3,'--k')

xlabel('q (‰)')
ylabel('p (‰)')

hold on
q1=-0.5;
p1=0.65;
sz=150;
scatter(q1,p1,sz,'.k')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This section of code plots model results for specified values of p and q.
% Values of q1 and p1 (lines 262-263) can be changed depending on
% the zone of overlap generated in Figure #1.

sz1 = 40;

eb1 = 2*d_4_in - 2*q1 ;
fc1 = (5*d_5_in - 2*eb1 - 5*q1 - 5*p1)/eb1 ;
fd1 = (7*d_7_in - 7*q1 - 21*p1 - 2*eb1 - 2*fc1*eb1)/eb1 ;
fe1 = (8*S8_HS_offset - 2*eb1 - 2*fc1*eb1 - 2*fd1*eb1) / eb1 ;

post_d_4_a = q1 ;
post_d_4_b = q1 + eb1 ;

post_d_5_a = q1 + p1 ;
post_d_5_b = q1 + p1 + eb1 ;
post_d_5_c = q1 + p1 + fc1*eb1 ;

post_d_6_a = q1 + 2*p1;
post_d_6_b = q1 + 2*p1 + eb1;
post_d_6_c = q1 + 2*p1 + fc1*eb1;

post_d_7_a = q1 + 3*p1 ;
post_d_7_b = q1 + 3*p1 + eb1 ;
post_d_7_c = q1 + 3*p1 + fc1*eb1 ;
post_d_7_d = q1 + 3*p1 + fd1*eb1 ;

post_d_8_a = q1 + 4*p1 ;
post_d_8_b = q1 + 4*p1 + eb1 ;
post_d_8_c = q1 + 4*p1 + fc1*eb1 ;
post_d_8_d = q1 + 4*p1 + fd1*eb1 ;

post_d_9_a = q1 + 5*p1 ;
post_d_9_b = q1 + 5*p1 + eb1 ;
post_d_9_c = q1 + 5*p1 + fc1*eb1 ;
post_d_9_d = q1 + 5*p1 + fd1*eb1 ;
post_d_9_e = q1 + 5*p1 + fe1*eb1 ;




figure(2)
plot([1 2],[post_d_4_a post_d_4_b],'k')
hold on
scatter([1 2],[post_d_4_a post_d_4_b],sz1,'k','filled')
hold on
ylabel('{\delta}^{34}S offset from system average (‰)')
xlabel('Position in polysulfide')
xlim([0 6.5])
ylim([-2 12.5])
xticks([1 2 3 4 5])
xticklabels({'a','b','c','d','e'})

hold on
plot([1 2 3],[post_d_5_a post_d_5_b post_d_5_c],'k')
hold on
scatter([1 2 3],[post_d_5_a post_d_5_b post_d_5_c],sz1,'k','filled')
hold on
scatter([1 2],[post_d_5_a post_d_5_b],sz1,'r','filled')
hold on

d_pyr1 = (post_d_5_a + post_d_5_b)/2;
yline(d_pyr1)

plot([1 2 3],[post_d_6_a post_d_6_b post_d_6_c],'k')
hold on
scatter([1 2 3],[post_d_6_a post_d_6_b post_d_6_c],sz1,'k','filled')
hold on

plot([1 2 3 4],[post_d_7_a post_d_7_b post_d_7_c post_d_7_d],'k')
hold on
scatter([1 2 3 4],[post_d_7_a post_d_7_b post_d_7_c post_d_7_d],sz1,'k','filled')
hold on

plot([1 2 3 4],[post_d_8_a post_d_8_b post_d_8_c post_d_8_d],'k')
hold on
scatter([1 2 3 4],[post_d_8_a post_d_8_b post_d_8_c post_d_8_d],sz1,'k','filled')
hold on

plot([1 2 3 4 5],[post_d_9_a post_d_9_b post_d_9_c post_d_9_d post_d_9_e],'k')
hold on

scatter([1 2 3 4 5],[post_d_9_a post_d_9_b post_d_9_c post_d_9_d post_d_9_e],sz1,'k','filled')
hold on

d_s81 = (post_d_9_a + 2*post_d_9_b + 2*post_d_9_c + 2*post_d_9_d + post_d_9_e)/8;
yline(d_s81)
