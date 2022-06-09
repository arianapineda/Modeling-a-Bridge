% Ariana Pineda, CAAM 210, SPRING 2022, Bridge Project
% bridge_1_driver.m
% this script models a simple truss bridge of arbitrary length using a 
% dynamically created adjacency matrix and plots deformed bridge with 
% weight of several cars
% Last modified: March 23, 2022

function bridge_1_driver()
%TASK 1

%plots spy of adjacency matrix
figure
plot_spy(1)
figure
plot_spy(2)
figure
plot_spy(3)

%plots the bridge with  0 load
figure
[~ , xc, yc, ~] = build_basic_bridge(1);
plot_bridge(1,xc,yc);
figure
[~ , xc, yc, ~] = build_basic_bridge(2);
plot_bridge(2,xc,yc);
figure
[~ , xc, yc, ~] = build_basic_bridge(3);
plot_bridge(3,xc,yc);

% TASK 3

%plots deformed bridge with load of 0.01 units
figure
build_load_plot_basic_bridge(1, 0.01);
figure
build_load_plot_basic_bridge(2, 0.01);
figure
build_load_plot_basic_bridge(3, 0.01);

%plots deformed bridge with load of 0.05 units
figure
build_load_plot_basic_bridge(1, 0.05);
figure
build_load_plot_basic_bridge(2, 0.05);
figure
build_load_plot_basic_bridge(3, 0.05);

end


function plot_bridge(nos,xc,yc)
% This function plots a basic bridge
% Inputs: nos—number of sections, xc—x coordinates of nodes; yc—y
% coordinates of nodes

%plots bridge structure
hold on
line(xc',yc','Color','b')
title([num2str(nos),' Section Bridge with Zero Load']);

%removes axis ticks
set(gca,'xtick',[ ]);
set(gca,'ytick',[ ]);
hold off
hold on

%plots base structures
fill([0 -2 -2 0],[0 0 -1 -1],'k');
fill([nos+2 nos+4 nos+4 nos+2],[0 0 -1 -1],'k');

%set plot bounds
xlim([-2, nos+4]);
ylim([-1 1.2]);
hold off
end

function plot_spy(nos)
% This function plots a spy plot of the adjacency matrix
% Inputs: nos—number of sections
[adj]=build_basic_bridge(nos);
spy(adj)

%labeling
title(nos + " Section Bridge Adjacency Matrix");
xlabel("nonzero matrix columns");
ylabel("nonzero matrix rows");
end

function build_load_plot_basic_bridge(nos, car_weight)
% this function plots the deformed bridge
% Inputs: nos—number of sections, car_weight—weight of the car

%solves for force in terms of applied weight
f = 5*nos+5;
n = 2*nos+2;
force = zeros(2*n, 1);
force(2,1) = -car_weight;

%fills in values of force vecor
for i=1:nos
    force(i*4+2,1) = -car_weight;
end

%get outputs of deform_basic_bridge function
[adj, xc, yc, len] = build_basic_bridge(nos);
[dx,dy,work,~,~]=deform_basic_bridge(nos,adj,xc,yc,len,force);

%remove axis ticks
set(gca,'xtick',[ ]);
set(gca,'ytick',[ ]);

%plots deformed bridge for each fiber
for i=1:f
    hold on
    plot(dx(i, :), dy(i, :), "b-");
    fill([0 -2 -2 0],[0 0 -1 -1],'k');
    fill([nos+2 nos+4 nos+4 nos+2],[0 0 -1 -1],'k');

end

%title
title(nos+" section bridge when car's weight " + car_weight + " units (work = " + work + ")");
hold off
end

function [adj, xc, yc, len] = build_basic_bridge(NoS)
% this function builds a bridge with 0 force applied
% Inputs: nos—number of sections
% Outputs: adj-adjacency matrix; xc-x coordinates of nodes; yc-y
% coordinates of nodes; len-fiber lengths

% calculate some helpful numbers
fibers = 5*NoS+5;
Nodes=2+2*NoS;
s=1/sqrt(2);
% initialize the return values
adj=zeros(fibers,Nodes*2);
xc = zeros(fibers, 2);
yc = zeros(fibers, 2);
len = ones (fibers, 1);
% build the left side of bridge
adj(1,1) = 1;
adj(2,[3,4]) = [s s];
xc(1, :) = [0 1];
xc(2, :) = [0 1];
yc(1, :) = [0 0];
yc(2, :) = [0 1];
len(2) = 1/s;

%build middle using a for loop
for n=0:NoS-1
    %sets starting row and col
    R=5*n;
    C=n*4+1;

    %calculate coordinates
    xc(R+3, :) = [n+1 n+1];
    yc(R+3, :) = [0 1];
    xc(R+4, :) = [n+1 n+2];
    yc(R+4, :) = [1 0];
    xc(R+5, :) = [n+1 n+2];
    yc(R+5, :) = [1 1];
    xc(R+6, :) = [n+1 n+2];
    yc(R+6, :) = [0 1];
    xc(R+7, :) = [n+1 n+2];
    yc(R+7, :) = [0 0];

    %fill in adjacency matrix
    adj(R+3,[C+0 C+1 C+2 C+3])=[0 -1 0 1];
    adj(R+4,[C+2 C+3 C+4 C+5])=[-1/s 1/s 1/s -1/s];
    adj(R+5,[C+2 C+3 C+6 C+7])=[-1 0 1 0];
    adj(R+6,[C+0 C+1 C+6 C+7])=[-1/s -1/s 1/s 1/s];
    adj(R+7,[C+0 C+1 C+4 C+5])=[-1 0 1 0];

end

% build the right side of bridge

%set starting row and col
start_col = (NoS*4)+1;
start_row=5*NoS+3;

%fill in adjacency matrix
adj(start_row,[start_col+0 start_col+1 start_col+2 start_col+3]) = [0 -1 0 1];
adj(start_row+1,[start_col+2 start_col+3]) = [1/s -1/s];
adj(start_row+2,[start_col+0 start_col+1]) = [1 0];

%fill in coordinates
xc(start_row, :) = [NoS+1 NoS+1];
yc(start_row, :) = [0 1];
xc(start_row+1, :) = [NoS+1 NoS+2];
yc(start_row+1, :) = [1 0];
xc(start_row+2, :) = [NoS+1 NoS+2];
yc(start_row+2, :) = [0 0];
end

function [dx, dy, work, X, Y] = deform_basic_bridge(nos,adj,xc,yc,len,F)
% this function builds a deformed bridge
% Inputs: nos—number of sections; adj-adjacency matrix; xc-x coordinates of
% nodes; yc-y coordinates of nodes; len-fiber lengths; F-force applied
% Outputs: dx-deformed x coord matrix; dy-deformed y coord matrix; work-work;
% X-x displacements vector; Y-y displacements vector

% calculate some helpful numbers
stiffness = adj'*diag(1./len)*adj;
displacements = stiffness\F;
work = displacements'*F;

% initialize the return values
X = displacements(1:2:end);
Y = displacements(2:2:end);
dx = zeros(size(xc));
dy = zeros(size(yc));

% deform the left side of bridge
dx(1,:) = xc(1,:) + [0 X(1)] ;
dx(2,:) = xc(2,:) + [0 X(2)] ;
dy(1,:) = yc(1,:) + [0 Y(1)] ;
dy(2,:) = yc(2,:) + [0 Y(2)];

% deform the middle of bridge
for n=0:nos-1

    %set starting row and col
    r=5*(n+1);
    c = 2*(n+1);

    %fill in values of dx and dy matrix
    dx(r,:)=xc(r-2,:)+[X(c-),X(c-2)];
    dy(r,:)=yc(r-2,:)+[Y(c-3) Y(c)];

    dx(r+1,:)=xc(r+1,:)+[X(c),X(c+1)];
    dy(r+1,:)=yc(r+1,:)+[Y(c) Y(c+1)];

    dx(r+2,:)=xc(r+2,:)+[X(c),X(c+2)];
    dy(r+2,:)=yc(r+2,:)+[Y(c) Y(c+2)];

    dx(r+3,:)=xc(r+3,:)+[X(c-1),X(c+2)];
    dy(r+3,:)=yc(r+3,:)+[Y(c-1) Y(c+2)];

    dx(r+4,:)=xc(r+4,:)+[X(c-1),X(c+1)];
    dy(r+4,:)=yc(r+4,:)+[Y(c-1) Y(c+1)];

end

%deform right side of bridge

%initialize variables
r=5*nos+3;
num_nodes=2+2*nos;

%fill in values of dx and dy matrix
dx(r,:)=xc(r,:)+[X(num_nodes-1),X(num_nodes)];
dy(r,:)=yc(r,:)+[Y(num_nodes-1),Y(num_nodes)];

dx(r+1,:)=xc(r+1,:)+[X(num_nodes),0];
dy(r+1,:)=yc(r+1,:)+[Y(num_nodes),0];

dx(r+2,:)=xc(r+2,:)+[X(num_nodes-1),0];
dy(r+2,:)=yc(r+2,:)+[Y(num_nodes-1),0];
end

% 1. My bridge should not reach below -0.2 on the y axis in order to be
% considered safe for a light-weight car . Therefore, 3 sections is the
% maximum length to fit this criteria.

% 2. My bridge should not reach below -0.2 on the y-axis in order to be
% considered safe for a car weighing 0.05 units. Therefore the maximum
% number of sections the bridge can have is one.

% 3. If there is more deformation in the bridge, the number of dots that
% appear on the SPY plot increases because there are less areas where the
% x or y coordinates of nodes relative to a fiber is 0. We would expect
% less fibers to be perfectly horizontal or vertical.

% 4. When nos increases, the bridge deforms more and saggs more. This is
% because there are more areas available for the force to push down on.
% Additionally, there are more dots that show up on the SPY.

