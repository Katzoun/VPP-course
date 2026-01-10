clear;clc;close all;
GW = createGridWorld(3,4);
GW.CurrentState = "[2,2]";
GW.TerminalStates = ["[1,4]"];
GW.ObstacleStates = ["[2,2]";"[2,4]"];

% u=1 - north, u=2 - south, u=3 - east, u=4 - west

%% north action
GW.T(state2idx(GW,"[1,1]"),[state2idx(GW,"[1,1]"),state2idx(GW,"[1,2]")], 1) = [0.9,0.1]; 
GW.T(state2idx(GW,"[1,2]"),[state2idx(GW,"[1,2]"),state2idx(GW,"[1,1]"),state2idx(GW,"[1,3]")],1) = [0.8,.1,.1]; 
GW.T(state2idx(GW,"[1,3]"),[state2idx(GW,"[1,3]"),state2idx(GW,"[1,2]"),state2idx(GW,"[1,4]")],1) = [0.8,.1,.1]; 
GW.T(state2idx(GW,"[2,1]"),[state2idx(GW,"[1,1]"),state2idx(GW,"[2,1]")],1) = [0.8,.2]; 
GW.T(state2idx(GW,"[2,3]"),[state2idx(GW,"[1,3]"),state2idx(GW,"[2,3]")],1) = [0.8,.2];
GW.T(state2idx(GW,"[2,4]"),[state2idx(GW,"[1,4]"),state2idx(GW,"[2,4]")],1) = [0,1];
GW.T(state2idx(GW,"[3,1]"),[state2idx(GW,"[3,1]"),state2idx(GW,"[3,2]"),state2idx(GW,"[2,1]")],1) = [0.1,0.1,0.8]; 
GW.T(state2idx(GW,"[3,2]"),[state2idx(GW,"[3,2]"),state2idx(GW,"[3,1]"),state2idx(GW,"[3,3]"),state2idx(GW,"[2,2]")],1) = [0.8,.1,.1,0]; 
GW.T(state2idx(GW,"[3,3]"),[state2idx(GW,"[2,3]"),state2idx(GW,"[3,2]"),state2idx(GW,"[3,4]")],1) = [0.8,.1,.1]; 
GW.T(state2idx(GW,"[3,4]"),[state2idx(GW,"[3,4]"),state2idx(GW,"[3,3]"),state2idx(GW,"[2,4]")],1) = [0.9,0.1,0]; 

%% south action
GW.T(state2idx(GW,"[1,1]"),[state2idx(GW,"[1,1]"),state2idx(GW,"[2,1]"),state2idx(GW,"[1,2]")],2) = [0.1,0.8,0.1]; 
GW.T(state2idx(GW,"[1,2]"),[state2idx(GW,"[1,2]"),state2idx(GW,"[1,1]"),state2idx(GW,"[1,3]"),state2idx(GW,"[2,2]")],2) = [0.8,0.1,0.1,0];
GW.T(state2idx(GW,"[1,3]"),[state2idx(GW,"[2,3]"),state2idx(GW,"[1,2]"),state2idx(GW,"[1,4]")],2) = [0.8,0.1,0.1]; 
GW.T(state2idx(GW,"[2,1]"),[state2idx(GW,"[3,1]"),state2idx(GW,"[2,1]")],2) = [0.8,0.2]; 
GW.T(state2idx(GW,"[2,3]"),[state2idx(GW,"[3,3]"),state2idx(GW,"[2,3]")],2) = [0.8,.2]; 
GW.T(state2idx(GW,"[3,1]"),[state2idx(GW,"[3,1]"),state2idx(GW,"[3,2]")],2) = [0.9,0.1]; 
GW.T(state2idx(GW,"[3,2]"),[state2idx(GW,"[3,2]"),state2idx(GW,"[3,1]"),state2idx(GW,"[3,3]")],2) = [0.8,0.1,0.1]; 
GW.T(state2idx(GW,"[3,3]"),[state2idx(GW,"[3,3]"),state2idx(GW,"[3,2]"),state2idx(GW,"[3,4]")],2) = [0.8,0.1,0.1]; 
GW.T(state2idx(GW,"[3,4]"),[state2idx(GW,"[3,4]"),state2idx(GW,"[3,3]")],2) = [0.9,0.1]; 

%% action east
GW.T(state2idx(GW,"[1,1]"),[state2idx(GW,"[1,1]"),state2idx(GW,"[1,2]"),state2idx(GW,"[2,1]")],3) = [0.1,0.8,0.1]; 
GW.T(state2idx(GW,"[1,2]"),[state2idx(GW,"[1,1]"),state2idx(GW,"[1,3]"),state2idx(GW,"[2,2]")],3) = [0.2,0.8,0]; 
GW.T(state2idx(GW,"[1,3]"),[state2idx(GW,"[1,3]"),state2idx(GW,"[1,4]"),state2idx(GW,"[2,3]")],3) = [0.1,0.8,0.1]; 
GW.T(state2idx(GW,"[2,1]"),[state2idx(GW,"[2,1]"),state2idx(GW,"[1,1]"),state2idx(GW,"[3,1]"),state2idx(GW,"[2,2]")],3) = [0.8,0.1,0.1,0]; 
GW.T(state2idx(GW,"[2,3]"),[state2idx(GW,"[2,3]"),state2idx(GW,"[1,3]"),state2idx(GW,"[3,3]"),state2idx(GW,"[2,4]")],3) = [0.8,0.1,0.1,0]; 
GW.T(state2idx(GW,"[3,1]"),[state2idx(GW,"[3,1]"),state2idx(GW,"[3,2]"),state2idx(GW,"[2,1]")],3) = [0.1,0.8,0.1]; 
GW.T(state2idx(GW,"[3,2]"),[state2idx(GW,"[3,3]"),state2idx(GW,"[3,2]")],3) = [0.8,0.2]; 
GW.T(state2idx(GW,"[3,3]"),[state2idx(GW,"[3,3]"),state2idx(GW,"[3,4]"),state2idx(GW,"[2,3]")],3) = [0.1,0.8,0.1]; 
GW.T(state2idx(GW,"[3,4]"),[state2idx(GW,"[3,4]")],3) = [1]; 

%% action west
GW.T(state2idx(GW,"[1,1]"),[state2idx(GW,"[1,1]"),state2idx(GW,"[2,1]")],4) = [0.9,0.1]; 
GW.T(state2idx(GW,"[1,2]"),[state2idx(GW,"[1,1]"),state2idx(GW,"[1,2]")],4) = [0.8,0.2]; 
GW.T(state2idx(GW,"[1,3]"),[state2idx(GW,"[1,2]"),state2idx(GW,"[1,3]"),state2idx(GW,"[2,3]")],4) = [0.8,0.1,0.1]; 
GW.T(state2idx(GW,"[2,1]"),[state2idx(GW,"[2,1]"),state2idx(GW,"[1,1]"),state2idx(GW,"[3,1]")],4) = [0.8,0.1,0.1]; 
GW.T(state2idx(GW,"[2,3]"),[state2idx(GW,"[2,3]"),state2idx(GW,"[1,3]"),state2idx(GW,"[3,3]"),state2idx(GW,"[2,2]")],4) = [0.8,0.1,0.1,0]; 
GW.T(state2idx(GW,"[3,1]"),[state2idx(GW,"[3,1]"),state2idx(GW,"[2,1]")],4) = [0.9,0.1]; 
GW.T(state2idx(GW,"[3,2]"),[state2idx(GW,"[3,1]"),state2idx(GW,"[3,2]")],4) = [0.8,0.2]; 
GW.T(state2idx(GW,"[3,3]"),[state2idx(GW,"[3,2]"),state2idx(GW,"[3,1]"),state2idx(GW,"[2,3]")],4) = [0.8,0.1,0.1]; 
GW.T(state2idx(GW,"[3,4]"),[state2idx(GW,"[3,3]"),state2idx(GW,"[3,4]")],4) = [0.8,0.2]; 
vec = zeros(1,12); vec(5) = 1;
GW.T(5,:,1) = vec; GW.T(5,:,2) = vec; GW.T(5,:,3) = vec; GW.T(5,:,4) = vec;

env = rlMDPEnv(GW);
plot(env);hold on ;
text(4,-1,'0','FontSize',15,'HorizontalAlignment','center');

%% value iter
Js = 4*ones(12,1); 
Js(11) = 0; Js(10) = 0; mu = ones(12,1); epsilon = 0.0000001;
J_new = zeros(12,1);
for i=1:2000
    J_new = zeros(12,1);
    for j=1:12
        if j == 10 || j == 11 || j == 5
            J_new(j) = 0;
        else
        [minval,minpos] = min([1+GW.T(j,:,1)*Js(:,i),...
                        1+GW.T(j,:,2)*Js(:,i),...
                        1+GW.T(j,:,3)*Js(:,i),...
                        1+GW.T(j,:,4)*Js(:,i)]);
        J_new(j) = minval; mu(j) = minpos;
        end
    end
    c_l = min(J_new-Js(:,i)); c_u = max(J_new-Js(:,i));
    if norm(c_l-c_u) < epsilon
        break;
    end
    Js = [Js,J_new];
end
J_new

%% policy iter
% mu = ones(12,1); mu_new = zeros(12,1); iter = 0;
% while 1
%     iter = iter  + 1;
%     for i=1:12
%         if i == 5 || i == 10 || i == 11
%             P(i,i) = 0;
%             g(i,1) = 0;
%         else
%             P(i,:) = GW.T(i,:,mu(i));
%             g(i,1) = 1;
%         end
%     end
%     J = (eye(12)-P)\g;
%     for i=1:12
%         if i == 5 || i == 10 || i == 11
%             mu_new(i) = 1;
%         else
%             [minval,minpos] = min([1+GW.T(i,:,1)*J,...
%                          1+GW.T(i,:,2)*J,...
%                          1+GW.T(i,:,3)*J,...
%                          1+GW.T(i,:,4)*J]);
%             mu_new(i) = minpos;
%         end
%     end
%     if all(mu_new == mu)
%         break;
%     else
%         mu = mu_new;
%     end
% end
% J_new = J;

%% print solution
strs = ['↑','↓','→','←'];
text(1,-1.4,strs(mu(1)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(1,-1,num2str(J_new(1)),'FontSize',15,'HorizontalAlignment','center');
text(1,-2.4,strs(mu(2)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(1,-2,num2str(J_new(2)),'FontSize',15,'HorizontalAlignment','center');
text(1,-3.4,strs(mu(3)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(1,-3,num2str(J_new(3)),'FontSize',15,'HorizontalAlignment','center');
text(2,-1.4,strs(mu(4)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(2,-1,num2str(J_new(4)),'FontSize',15,'HorizontalAlignment','center');
text(2,-3.4,strs(mu(6)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(2,-3,num2str(J_new(6)),'FontSize',15,'HorizontalAlignment','center');
text(3,-1.4,strs(mu(7)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(3,-1,num2str(J_new(7)),'FontSize',15,'HorizontalAlignment','center');
text(3,-2.4,strs(mu(8)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(3,-2,num2str(J_new(8)),'FontSize',15,'HorizontalAlignment','center');
text(3,-3.4,strs(mu(9)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(3,-3,num2str(J_new(9)),'FontSize',15,'HorizontalAlignment','center');
text(4,-3.4,strs(mu(12)),'FontSize',20,'Color','k','HorizontalAlignment','center','VerticalAlignment','middle');
text(4,-3,num2str(J_new(12)),'FontSize',15,'HorizontalAlignment','center');

% figure
% plot(Js(1,:))