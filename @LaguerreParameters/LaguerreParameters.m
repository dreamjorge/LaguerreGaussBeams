classdef LaguerreParameters < GaussianParameters & handle & matlab.mixin.Copyable
% LaguerreParameters class gives Laguerre parameters of Laguerre Gaussian Beam
% recives zCoordinate,InitialWaist,Wavelength,l,p as input and
% gives a object with next properties:
%
% -l
% -p
% -LaguerreWaist
% -PhiPhase
% -zCoordinate
% -InitialWaist
% -Wavelength
% -RayleighDistance
% -k
% -Waist
% -Radius
% -GouyPhase
% -Amplitude
% -DivergenceAngle
% 

% Example: Parameters = LaguerreParameters(zCoordinate...
%                                         ,InitialWaist...
%                                         ,Wavelength)
  
  properties
    l % angular number
    p % radial number
  end
  
  properties (Dependent)
    laguerreWaist  % Waist of Laguerre Gauss Beam
    phiPhase       % Phase of Laguerre Gauss Beam
  end
  
  methods(Static)
      
    function waistL = getwaist(zCoordinate,InitialWaist,RayleighDistance,p,l)
    %%Function for estimate wast of Laguerre Gauss Beam
      waistL = sqrt((2*p+l+1))...
             .*(InitialWaist).*sqrt((zCoordinate./RayleighDistance).^2+1);  
    end
     

    
  end
  
  
  methods
    
    function PhiPhase      = get.phiPhase(obj) 
      % Function for estimate phase of Laguerre Gaussian Beam
      PhiPhase = (abs(obj.l)+2*(obj.p)-1).*obj.GouyPhase; 
    end

    function LaguerreWaist = get.laguerreWaist(obj)
      % Function for estimate phase of Laguerre Gaussian Beam
      LaguerreWaist = LaguerreParameters.getwaist(obj.zCoordinate,...
                                                  obj.initialWaist,obj.RayleighDistance,...
                                                  obj.l,...
                                                  obj.p);
    end

    
    function Parameters = LaguerreParameters(zCoordinate, ...
                                             initialWaist, ...
                                             wavelength, ...
                                             l, ...
                                             p, ...
                                             units)
      %% LaguerreParameters Object
      % generation of object with arguments inputs
      
%       arguments
%         zCoordinate   (1,:) double {mustBeNumeric,mustBeReal}
%         initialWaist  (1,1) double {mustBeNumeric,mustBeReal,mustBePositive}
%         wavelength    (1,1) double {mustBeNumeric,mustBeReal,mustBePositive}
%         l             (1,1) double {mustBeNumeric,mustBeReal,mustBeInteger}
%         p             (1,1) double {mustBeNumeric,mustBeReal,mustBeInteger}
%         units         (1,1) string
%       end
      %Call Gaussian Parameters
      Parameters@GaussianParameters(zCoordinate,...
                                    initialWaist,...
                                    wavelength,...
                                    units);

      % Add l,p parameters of input to object
      Parameters.l = l;
      Parameters.p = p;


    end
    
    
    function [] = plotParameters(obj)
      %% Function plots parameters of Laguerre Gaussian Beam
      % Input:
      %  
      p1            = plot(obj.zCoordinate ,obj.waist,'Color','red');
      hold on
      p2            = plot(obj.zCoordinate,-obj.waist,'Color','red');
      p3            = plot(obj.zCoordinate, obj.zCoordinate*tan(obj.divergenceAngle),'Color','blue');
      p4            = plot(obj.zCoordinate,-obj.zCoordinate*tan(obj.divergenceAngle),'Color','blue');
      p5            = plot(obj.zCoordinate,-obj.laguerreWaist,'Color','green');
      p6            = plot(obj.zCoordinate, obj.laguerreWaist,'Color','green');
      xlabel(['Distance of Propagation on ',obj.units])
      ylabel(['x on ',obj.units])
      title('Parameters of Laguerre Gaussian Beam')
      legend([p1,p5,p3],{'Waist of Gaussian Beam'...
                        ,'Waist of Laguerre Beam'...
                        ,['Angle of Divergence = ',num2str(rad2deg(obj.divergenceAngle)),'°']})
      hold off
      
 
      
      end

   end
end