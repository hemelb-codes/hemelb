package uk.ac.ucl.chem.ccs.vizclient;

import java.io.Serializable;

public class SteeringData implements Serializable {
	
	private float ctr_x;
	private float ctr_y;
	private float ctr_z;
	private float longitude;
	private float latitude;
	private float zoom_factor;
	private float vis_brightness;
	private float velocity_max;
	private float stress_max;
	private int vis_mouse_x;
	private int vis_mouse_y;
	private int lbm_terminate;
	private float pressure_min;		
	private float pressure_max;		
	private float vis_glyph_length;	
	private int pixels_x;			
	private int pixels_y;			
	private int vis_mode;
	private float vis_streaklines_per_pulsatile_period;
	private float vis_streakline_length;
	
	
//	public SteeringData(float ctr_x, float ctr_y, float ctr_z, float latitude,
//			float longitude, float stress_max, float velocity_max,
//			float vis_brightness, float zoom_factor, int vis_mouse_x, int vis_mouse_y, int kill) {
//		super();
//		this.ctr_x = ctr_x;
//		this.ctr_y = ctr_y;
//		this.ctr_z = ctr_z;
//		this.latitude = latitude;
//		this.longitude = longitude;
//		this.stress_max = stress_max;
//		this.velocity_max = velocity_max;
//		this.vis_brightness = vis_brightness;
//		this.zoom_factor = zoom_factor;
//		
//		
//		
//		this.vis_mouse_x = vis_mouse_x;
//		this.vis_mouse_y = vis_mouse_y;
//		this.lbm_terminate = kill;
//	}
	
	public SteeringData() {
		super();
		this.ctr_x = 0.0f;
		this.ctr_y = 0.0f;
		this.ctr_z = 0.0f;
		this.latitude = 45.0f;
		this.longitude = 45.0f;
		this.zoom_factor = 1.0f;
		this.vis_brightness = 0.03f;
		this.velocity_max = 0.1f;
		this.stress_max = 0.1f;
		this.vis_mouse_x = -1;
		this.vis_mouse_y = -1;
		this.lbm_terminate = 0;
		this.pressure_min = 80.0f;		
		this.pressure_max = 120.0f;		
		this.vis_glyph_length = 1.0f;	
		this.pixels_x = 512;			
		this.pixels_y = 512;			
		this.vis_mode = 0;
		this.vis_streaklines_per_pulsatile_period = 5.0f;
		this.vis_streakline_length = 100.0f;	
	}
	
	public float getPressure_min() {
		return pressure_min;
	}


	public void setPressure_min(float pressure_min) {
		this.pressure_min = pressure_min;
	}

	public float getPressure_max() {
		return pressure_max;
	}



	public void setPressure_max(float pressure_max) {
		this.pressure_max = pressure_max;
	}



	public float getVis_glyph_length() {
		return vis_glyph_length;
	}



	public void setVis_glyph_length(float vis_glyph_length) {
		this.vis_glyph_length = vis_glyph_length;
	}



	public int getPixels_x() {
		return pixels_x;
	}



	public void setPixels_x(int pixels_x) {
		this.pixels_x = pixels_x;
	}



	public int getPixels_y() {
		return pixels_y;
	}



	public void setPixels_y(int pixels_y) {
		this.pixels_y = pixels_y;
	}



	public int getVis_mode() {
		return vis_mode;
	}



	public void setVis_mode(int vis_mode) {
		this.vis_mode = vis_mode;
	}



	public float getVis_streaklines_per_pulsatile_period() {
		return vis_streaklines_per_pulsatile_period;
	}



	public void setVis_streaklines_per_pulsatile_period(
			float vis_streaklines_per_pulsatile_period) {
		this.vis_streaklines_per_pulsatile_period = vis_streaklines_per_pulsatile_period;
	}



	public float getVis_streakline_length() {
		return vis_streakline_length;
	}



	public void setVis_streakline_length(float vis_streakline_length) {
		this.vis_streakline_length = vis_streakline_length;
	}



	public int getVis_mouse_x() {
		return vis_mouse_x;
	}

	public void setVis_mouse_x(int vis_mouse_x) {
		this.vis_mouse_x = vis_mouse_x;
	}

	public int getVis_mouse_y() {
		return vis_mouse_y;
	}

	public void setVis_mouse_y(int vis_mouse_y) {
		this.vis_mouse_y = vis_mouse_y;
	}

	public int getKill() {
		return lbm_terminate;
	}

	public void setKill(int kill) {
		this.lbm_terminate = kill;
	}

	public float getCtr_x() {
		return ctr_x;
	}

	public void setCtr_x(float ctr_x) {
		this.ctr_x = ctr_x;
	}

	public float getCtr_y() {
		return ctr_y;
	}

	public void setCtr_y(float ctr_y) {
		this.ctr_y = ctr_y;
	}

	public float getCtr_z() {
		return ctr_z;
	}

	public void setCtr_z(float ctr_z) {
		this.ctr_z = ctr_z;
	}

	public float getLongitude() {
		return longitude;
	}

	public void setLongitude(float longitude) {
		this.longitude = longitude;
	}

	public void updateLongitude(float newLongitude) {
		this.longitude = (longitude + newLongitude)%360;
	}
	
	public float getLatitude() {
		return latitude;
	}

	public void setLatitude(float latitude) {
		this.latitude = latitude;
	}
	
	public void updateLatitude(float newLatitude) {
				
		this.latitude = (latitude + newLatitude)%360;
	}

	public void updateLongitudeLatitude(double dx, double dy) {
		
		this.latitude += dy; 
		
	//	if( (int) Math.floor( Math.abs(this.latitude) + 90.0 ) % 2 == 0 ) {
	//		this.longitude -= dx;
	//	} else {
			this.longitude += dx;
		//}

			this.latitude = this.latitude%360;
			this.longitude = this.longitude%360;			
			
		//System.out.println("Latitude " + latitude + " Longitude " + longitude);
		

	}
	
	
	public float getZoom_factor() {
		return zoom_factor;
	}

	public void setZoom_factor(float zoom) {
		//check bounds
		if (zoom > 40.0f) {
			this.zoom_factor = 40.0f;
		} else if (zoom < 0.1f) {
			this.zoom_factor = 0.1f;
		} else {
			this.zoom_factor = zoom;
		}	
	}

	public float getVis_brightness() {
		return vis_brightness;
	}

	public void setVis_brightness(float vis_brightness) {
		this.vis_brightness = vis_brightness;
	}

	public float getVelocity_max() {
		return velocity_max;
	}

	public void setVelocity_max(float velocity_max) {
		this.velocity_max = velocity_max;
	}

	public float getStress_max() {
		return stress_max;
	}

	public void setStress_max(float stress_max) {
		this.stress_max = stress_max;
	}
	
	
}
