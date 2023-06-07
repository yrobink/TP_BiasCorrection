
is_leap = function(years)##{{{
{
	return( ( (years %% 4 == 0) & !(years %%100 == 0) ) | (years %% 400 == 0) )
}
##}}}

inverse_rel_dist_earth_Sun = function( time )##{{{
{
	days = as.POSIXlt( time , format = "%Y-%m-%d" )$yday
	day_in_year = 365 + is_leap( as.integer(format.Date( time , "%Y" )) )
	d_r = 1 + 0.033 * base::cos( 2 * pi * days / day_in_year )
	return(d_r)
}
##}}}

solar_declination = function( time )##{{{
{
	days = as.POSIXlt( time , format = "%Y-%m-%d" )$yday
	day_in_year = 365 + is_leap( as.integer(format.Date( time , "%Y" )) )
	delt = 0.409 * base::sin( 2 * pi * days / day_in_year - 1.39 )
	
	return(delt)
}
##}}}

sunset_hour_angle = function( lat , delt )##{{{
{
	vec = tan(delt)
	mat = tan(lat * pi / 180)
	arg = vec * mat
	w_s = acos(-arg)
	return(w_s)
}
##}}}


atmospheric_pressure = function( height )##{{{
{
	#' Atmospheric pressure (P)
	#' ========================
	#' 
	#' Rapport FAO page 31 Equation (7)
	#' Units : 
	#' - P: [kPa]
	#' - z: [m] (above sea level)
	sea_ps = 101.3
	ps     = sea_ps * ( (293 - 0.0065 * height) / 293 )^5.26
	return(ps)
}
##}}}


saturation_vapour_airtemp = function(tasm)##{{{
{
	#' Saturation vapour pressure at the air temperature 
	#' =================================================
	#' 
	#' Rapport FAO page 36 Equation (11)
	#' Units :
	#' - e_o: [kPa]
	#' - T: [K]
	C = tasm - 273.15
	swv = 0.6108 * base::exp( (17.27 * C) / (C + 237.3) )
	return(swv)
}
##}}}

saturation_vapour = function( tasmax , tasmin )##{{{
{
	#' Mean saturation vapour pressure 
	#' ===============================
	#' 
	#' Rapport FAO page 36 Equation (12)
	#' Units :
	#' - e_s: [kPa]
	#' - tasmin, tasmax: [K]
	val_max = saturation_vapour_airtemp(tasmax)
	val_min = saturation_vapour_airtemp(tasmin)
	e_s = ( val_max + val_min ) / 2.
	return(e_s)
}
##}}}

saturation_vapour_pressure_tetens = function( tas )##{{{
{
	#' Saturation vapour pressure of water (Tetens) 
	#' ============================================
	#' 
	#' Wikipedia : https://en.wikipedia.org/wiki/Tetens_equation
	#' Units :
	#' - e_w: [kPa]
	#' - tas: [K]
	
	C   = tas - 273.15
	
	ewp = 0.6108 * base::exp(  17.27 * C / ( C + 237.3 ) )
	ewn = 0.6108 * base::exp( 21.875 * C / ( C + 265.5 ) )
	
	e_w = ewn * (tas < 273.15) + ewp * !(tas < 273.15)
	
	return(e_w)
}
##}}}

actual_saturation_h_esp = function( huss , ps )#{{{
{
	#' Actual vapour pressure (from specific humidity)
	#' ===============================================
	#' 
	#' https://fr.wikipedia.org/wiki/Humidit%C3%A9_sp%C3%A9cifique
	#' Units: 
	#' - e_part: [kPa]
	#' - altitud: [m] about sea level
	#' - h_e : [kg/kg]
	
	e_part = ps * huss / ( 0.622 + huss * 0.378 )
	return(e_part)
}
##}}}


extrater_radiation_daily = function( d_r , w , lat , delt )##{{{
{
	#' Extraterrestrial radiation for faily periods (Ra)
	#' =================================================
	#' 
	#' Rapport FAO page 46 Equation (21)
	#' Units : 
	#' - r_a: [MJ m^-2 day^-1]
	#' - d_r: [*]
	#' - w: [*]
	#' - lat: [rad] [*]
	#' - delt: [*]
	
	solar_constant = 0.0820
	ten1 = base::cos(delt) * base::cos(lat * pi / 180)
	ten2 = base::sin(delt) * base::sin(lat * pi / 180)
	ten3 = w * ten2 + ten1 * base::sin(w)
	r_a  = 24 * 60 * solar_constant * d_r * ten3 / pi
	return(r_a)
}
##}}}

clear_sky_solar_radiation = function( r_a , height )##{{{
{
	#' Clear-sky solar radiation 
	#' =========================
	#' 
	#' Rapport FAO page 51 Equation (37)
	#' Units : 
	#' - r_so : [MJ m^-2 day^-1]
	#' - z : [m] (above sea level)
	
	r_so = ( 0.75 + 0.00002 * height ) * r_a
	return(r_so)
}
##}}}


solar_radiation_derived = function( r_a , tasmax , tasmin , coef )##{{{
{
	#' Solar Radiation data derived from air temperature differences
	#' =============================================================
	#' 
	#' Rapport FAO page 60 Equation (50)
	#' Units : 
	#' - r_s: [MJ.m^-2.day^-1]
	#' 
	#' k_RS : Adjustment coefficient is 0.16 for interior locations and 0.19 for coastal location
	
	k_RS = coef
	r_s  = k_RS * r_a * (tasmax - tasmin)^0.5
	return(r_s)
}
##}}}

trans_unit_rayon = function(r_s)##{{{
{
	#' Transformation of radiation units
	#' =================================
	#' 
	#' Units : 
	#' - r_s input: [W.m^-2]
	#' - r_s output: [MJ.m^-2.day^-1]
	
	return(r_s * 24 * 60 * 60 / (1000 * 1000))
}
##}}}

net_outgoing_longwave_radiation = function( tasmax , tasmin , e_a , r_s , r_so )##{{{
{
	#' Net outgoing/longwave radiation
	#' ===============================
	#' 
	#' Rapport FAO page 52 Equation (39)
	#' Units : 
	#' - r_nl: [MJ m^-2 day^-1]
	#' - t: [K]
	#' r_s/r_so must be limited to <= 1.0
	
	div = r_s / r_so
	divp = div > 1
	
	if( base::any(divp) )
	{
		div[divp] = 1
	}
	
	r_nl = 4.903e-9 * ( ( tasmax^4+tasmin^4 ) / 2 ) * ( 0.34 - 0.139 * e_a^0.5 ) * ( 1.35 * div - 0.35 )
	
	return(r_nl)
}
##}}}

net_shortwave_radiation = function( r_s , albedo = 0.23 )##{{{
{
	#' Net solar/shortwave radiation
	#' =============================
	#' 
	#' Rapport FAO page 51 Equation (38)
	#' Units : 
	#' - r_ns: [MJ.m^-2.day^-1]
	#' - r_s: [MJ.m^-2.day^-1]
	
	r_ns = ( 1 - albedo ) * r_s
	return(r_ns)
}
##}}}


slope_saturation = function(tas)##{{{
{
	#' Slope of saturation vapour pressure curve (Delta)
	#' =================================================
	#' 
	#' Rapport FAO page 37 Equation (13)
	#' Units :
	#' - delta : [kPa°C^-1]
	#' - tas   : [K]
	
	C = tas - 273.15
	delta = 4098 * 0.6108 * base::exp( 17.27 * C / ( C + 237.3 ) ) / ( ( C + 237.3 )^2 )
	return(delta)
}
##}}}

psychrometric_cte = function(p)##{{{
{
	#' Psychrometric constant 
	#' ======================
	#' 
	#' Rapport FAO page 32 Equation (8)
	#' Units :
	#' - gamma: [kPa°C^-1]
	#' - pression: [kPa]
	
	gamma = 0.000665*p
	return(gamma)
}
##}}}

windSpeed2m = function( sfcWind , height_origin = 10 )##{{{
{
	#' Speed at 2 meters above ground surface
	#' ======================================
	#' 
	#' Rapport FAO page 56 Equation (47)
	#' Units 
	#' - sfcWind: [m s^-1]
	#' - w2m u_z: [m s^-1]
	#' - height_origin: [m] (from ground)
	
	w2m  = sfcWind * 4.87 / base::log( 67.8 * height_origin - 5.42 )
	w2m0 = (w2m < 0)
	if( base::any(w2m0) )
		w2m[w2m0] = 0.1
	
	return(w2m)
}
##}}}


penman_monteith_equation = function( delta , rn , g , gamma , tas , w2m , e_s , e_a )##{{{
{
	#' Penman Monteith Equation  
	#' ========================
	#' 
	#' Rapport FAO page 24 Equation (6)
	#' Units : [mm day^-1]
	
	dif  = e_s - e_a
	dif0 = !(dif > 0)
	if( base::any(dif0) )
		dif[dif0] = 0
	
	num = 0.408 * delta * ( rn - g ) + gamma * ( 900.0 / tas ) * w2m * dif
	den = delta + gamma * ( 1. + 0.34 * w2m )
	eto = num / den
	eton = !(eto > 0)
	if( base::any(eton) )
		eto[eton] = 0
	
	return(eto)
}
##}}}

human2SI_units = function(X)## {{{
{
	for( v in list("tas","tasmin","tasmax") )
		X[,v] = X[,v] + 273.15
	
	X[,"prtot"] = X[,"prtot"] / 86400
	
	return(X)
}
##}}}


find_ETP = function( X , zone = "Lez" , albedo = 0.23 , hargreaves_coef = 0.175 , weight = TRUE )##{{{
{
	X = human2SI_units(X)
	
	## Select latitude + height
	if( weight )
	{
		if( zone == "Lez" )
		{
			lat    = 43.73388247067313
			height = 144.25713471676863
		}
		
		if( zone == "Tech" )
		{
			lat    = 42.44337289551194
			height = 985.587075941803
		}
	}
	else
	{
		if( zone == "Lez" )
		{
			lat    = 43.73983414358141
			height = 125.5
		}
		
		if( zone == "Tech" )
		{
			lat    = 42.45164713469047
			height = 1118.857142857143
		}
	}
	
	## Extract variables
	time    = rownames(X)
	tasmin  = X[,"tasmin"]
	tasmax  = X[,"tasmax"]
	huss    = X[,"huss"]
	sfcWind = X[,"sfcWind"]
	tas     = (tasmin + tasmax) / 2
	
	## Now the (long) computation
	d_r  = inverse_rel_dist_earth_Sun(time)
	delt = solar_declination( time )
	w    = sunset_hour_angle( lat , delt )
	
	ps     = atmospheric_pressure(height)
	e_s    = saturation_vapour( tasmax , tasmin )
	e_w    = saturation_vapour_pressure_tetens(tas)
	e_part = actual_saturation_h_esp( huss , ps )
	e_a    = e_s * e_part / e_w
	
	r_a  = extrater_radiation_daily( d_r , w , lat , delt )
	r_so = clear_sky_solar_radiation( r_a , height )
	r_s  = solar_radiation_derived( r_a , tasmax , tasmin , hargreaves_coef )
	r_nl = net_outgoing_longwave_radiation( tasmax , tasmin , e_a , r_s , r_so )
	r_ns = net_shortwave_radiation( r_s , albedo )
	rn   = r_ns - r_nl
	
	delta = slope_saturation( tas )
	g     = 0
	gamma = psychrometric_cte(ps)
	w2m   = windSpeed2m(sfcWind)
	wpe   = penman_monteith_equation( delta , rn , g , gamma , tas , w2m , e_s , e_a )
	
	wpe   = matrix( wpe , ncol = 1 )
	rownames(wpe) = time
	colnames(wpe) = "evspsblpot"
	
	return(wpe)
}
##}}}


