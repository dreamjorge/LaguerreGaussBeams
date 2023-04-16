classdef GaussianParameters <  handle & matlab.mixin.Copyable
% This class gives Gaussian parameters of Gaussian Beam
% recives zCoordinate,InitialWaist,Wavelength as input and
% gives a object with next properties:
%
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

% Example: Parameters = GaussianParameters(zCoordinate,InitialWaist,Wavelength)

  properties(GetAccess =private)
    
  end

  properties
    %% Independent properties
    units            % physical units
    zCoordinate      % z-coordinate (or distance of propagation)
    initialWaist     % initial wasit of Gaussian Beam
    wavelength       % wavelength of Gaussian Beam
    
  end
  
  properties (Dependent)
    %% Properties dependent of Gaussian Beam
    RayleighDistance % RayleighDistance.
    k                % k-number.
    waist            % Waist of Gaussian Beam at zCoordinate.
    radius           % Radius of Gaussian Beam at zCoordinate.
    GouyPhase        % GouyPhace of Gaussian Beam at zCoordinate.
    amplitude        % Amplitude of Gaussian Beam at zCoordinate.
    divergenceAngle  % Angle of Divergence of Gaussian Beam
  end

  methods(Static)

    function waistZ = getwaist(zCoordinate,initialWaist,RayleighDistance) 
        
      waistZ = (initialWaist)*sqrt( (zCoordinate./RayleighDistance).^2 + 1);

    end
    
    function phiZ = getphase(zCoordinate,RayleighDistance)
        
      phiZ = atan(zCoordinate/RayleighDistance);
        
    end
    
    function radiusZ = getradius(zCoordinate,RayleighDistance)
        
      radiusZ = (zCoordinate).*(1+(RayleighDistance./zCoordinate).^2);

      if (find(isnan(radiusZ))~= 0)
          radiusZ(isnan(radiusZ))= inf;
      end
      
    end
                                             
  end
  methods

    function Parameters = GaussianParameters(zCoordinate, initialWaist, wavelength, units)
%       arguments
%         NameValueArgs.zCoordinate  (1,:) double {mustBeNumeric,mustBeReal}
%         NameValueArgs.initialWaist (1,1) double {mustBeNumeric,mustBeReal,mustBePositive}
%         NameValueArgs.wavelength   (1,1) double {mustBeNumeric,mustBeReal,mustBePositive}
%         NameValueArgs.units        (1,1) string
%       end
     %% Constructor of Gaussian Beam   
%       if nargin == 4 
        Parameters.zCoordinate  = zCoordinate;
        Parameters.initialWaist = initialWaist;
        Parameters.wavelength   = wavelength;
        Parameters.units        = units;
%       else
%          error('You need introduce zCoordinate, initialWaist, wavelength and units as inputs')
%       end
    end
    %% Error for dependent properties
    
    function [] = set.RayleighDistance(Parameters,~)
      fprintf('%s%d\n','RayleighDistance is: ',Parameters.RayleighDistance)
      error('You cannot set RayleighDistance property'); 
    end 
   
    function [] = set.k(Parameters,~)
      fprintf('%s%d\n','k is: ',Parameters.k)
      error('You cannot set k property'); 
    end 
    
    function [] = set.waist(Parameters,~)
      fprintf('%s%d\n','waist is: ',Parameters.Waist)
      error('You cannot set Waist property'); 
    end 
        
    function [] = set.radius(Parameters,~)
      fprintf('%s%d\n','radius is: ',Parameters.Radius)
      error('You cannot set Waist property'); 
    end 
       
    function [] = set.GouyPhase(Parameters,~)
      fprintf('%s%d\n','phase is: ',Parameters.PhiPhase)
      error('You cannot set PhiPhase property'); 
    end 
    
    function [] = set.amplitude(Parameters,~)
      fprintf('%s%d\n','amplitude is: ',Parameters.Amplitude)
      error('You cannot set Amplitude property'); 

    end 
    
    
    %% Get dependent properties
    
    function k = get.k(Parameters)
       k  = (2*pi)/Parameters.wavelength;    
    end
    
    function rayleighDistance = get.RayleighDistance(Parameters) 
      rayleighDistance  = ((Parameters.initialWaist).^2).*(Parameters.k)/2;    
    end
    
    function Waist = get.waist(Parameters)
      Waist = GaussianParameters.getwaist(Parameters.zCoordinate,...
                                          Parameters.initialWaist,...
                                          Parameters.RayleighDistance);
    end
    
    function Phase = get.GouyPhase(Parameters)
      Phase = GaussianParameters.getphase(Parameters.zCoordinate,...
                                          Parameters.RayleighDistance);
    end
    
    function Radius = get.radius(Parameters)
      Radius = GaussianParameters.getradius(Parameters.zCoordinate,...
                                            Parameters.RayleighDistance);
    end
    
    function Amplitude = get.amplitude(Parameters)
      Amplitude = 1./Parameters.waist;
    end
    
    function DivergenceAngle = get.divergenceAngle(Parameters)
      DivergenceAngle = atan((Parameters.initialWaist)./(Parameters.RayleighDistance));
    end

    function [fig] = plotParameters(Parameters)
    %% Function plots parameters of Gaussian Beam (Waist,DivergenceAngle)
    % Input:
    %  -GB Parameters as GaussianBeamParameters
      fig=figure;
      fig.Position=[320 263 1236 420];
      p1 = plot( Parameters.zCoordinate...
               , Parameters.waist...
               ,'Color','red' );
      hold on
      p2 = plot( Parameters.zCoordinate...
               ,-Parameters.waist...
               ,'Color','red' );  
      p3 = plot( Parameters.zCoordinate...
               , Parameters.zCoordinate*tan(Parameters.divergenceAngle)...
               ,'Color','blue');
      p4 = plot( Parameters.zCoordinate...
               ,-Parameters.zCoordinate*tan(Parameters.divergenceAngle)...
               ,'Color','blue');

      p5 = plot( Parameters.zCoordinate,...
                 Parameters.radius);      
      xlabel(['Distance of Propagation [',Parameters.units,']']);
      ylabel(['[',Parameters.units,']']);
      title('Parameters of Gaussian Beam');
      p1legend = 'Waist of Gaussian Beam';
      p3legend = ['Angle of Divergence = ',num2str(rad2deg(Parameters.divergenceAngle)),'°'];
      legend([p1,p3],{p1legend,p3legend});
      ylim([-250 250]);
    end

    
    
  end
  

end