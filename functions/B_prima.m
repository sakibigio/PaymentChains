% function that computes B' along a transition path
function [B_p,fval,flag] = B_prima(B,B_tilde,B_tilde_p,params)
    Bstar= params.Bstar;
    Bast = params.Bast;
    B_p_res=params.B_p_res;
  %  if B_tilde<=B_tilde_p
        if B<=Bstar(B_tilde)
            B_p=B; fval=0; flag=1;
        elseif B>=Bast(B_tilde_p)
            B_p=B; fval=0; flag=1;
        else
            if B>Bstar(B_tilde)
                [B_p,fval,flag]=fsolve(@(x) B_p_res(B,x,B_tilde,B_tilde_p),B_tilde_p-0.1);
                B_p=max(B_p,Bstar(B_tilde));
            else
                %if B>Bast(B_tilde,B_tilde_p)
                [B_p,fval,flag]=fsolve(@(x) B_p_res(B,x,B_tilde,B_tilde_p),Bstar(B_tilde)-0.1);
                B_p=max(B_p,Bstar(B_tilde));
            %elseif 
             %   [B_p,fval,flag]=fsolve(@(x) B_p_res(B,x,B_tilde,B_tilde_p),B_tilde_p-0.2);
                 if Bstar(B_tilde_p)<B
                   B_p=max(B_p,Bstar(B_tilde_p));
                 end
            end
        end
%     elseif B_tilde>B_tilde_p
%         if B>=Bstar(B_tilde)
%             B_p=B; fval=0; flag=1;
%         elseif B<=Bast(B_tilde_p)
%             B_p=B_tilde_p; fval=0; flag=1;
%         else
%             [B_p,fval,flag]=fsolve(@(x) B_p_res(B,x,B_tilde,B_tilde_p),Bstar(B_tilde+0.1));
%             B_p=max(B_p,Bstar(B_tilde));
%         end
%      elseif B_tilde==B_tilde_p
%          B_p=B; fval=0; flag=0;
   % else
    %    B_p=NaN; fval=NaN; flag=NaN;
    % end
end