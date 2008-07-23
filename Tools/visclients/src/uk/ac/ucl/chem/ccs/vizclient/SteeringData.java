package uk.ac.ucl.chem.ccs.vizclient;

public class SteeringData {
/*
	1. ctr_x - X pan - [0.0, 1.0] - for now just send 0.0 and no need to  
	steer this for now.
	2. ctr_y - Y pan - [0.0, 1.0] - for now just send 0.0 and no need to  
	steer this for now.
	3. ctr_z - Z pan - [0.0, 1.0] - for now just send 0.0 and no need to  
	steer this for now.
	4. longitude - conversion of left click and y mouse movement into  
	angle [-90.0,+90.0] degrees.
	5. latitude - conversion of left click and x mouse movement into angle  
	[-90.0,+90.0] degrees.
	6. zoom factor - [0.5 5.0] - Controllable using the mouse wheel - each  
	click on the wheel is a 0.1 increment
	7. vis_brightness - 0.03 - settable inside an editable field in a GUI  
	window
	8. velocity_max - 0.001 - settable inside an editable field in a GUI  
	window
	9. stress_max - 0.01 - settable inside an editable field in a GUI window
	*/
	
	private float ctr_x;
	private float ctr_y;
	private float ctr_z;
	private float longitude;
	private float latitude;
	private float zoom_factor;
	private float vis_brightness;
	private float velocity_max;
	private float stress_max;
	
	public SteeringData(float ctr_x, float ctr_y, float ctr_z, float latitude,
			float longitude, float stress_max, float velocity_max,
			float vis_brightness, float zoom_factor) {
		super();
		this.ctr_x = ctr_x;
		this.ctr_y = ctr_y;
		this.ctr_z = ctr_z;
		this.latitude = latitude;
		this.longitude = longitude;
		this.stress_max = stress_max;
		this.velocity_max = velocity_max;
		this.vis_brightness = vis_brightness;
		this.zoom_factor = zoom_factor;
	}
	
	public SteeringData() {
		super();
		this.ctr_x = 0.0f;
		this.ctr_y = 0.0f;
		this.ctr_z = 0.0f;
		this.latitude = 45.0f;
		this.longitude = 45.0f;
		this.zoom_factor = 1.0f;
		this.vis_brightness = 0.03f;
		this.velocity_max = 0.001f;
		this.stress_max = 0.001f;

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

	public float getZoom_factor() {
		return zoom_factor;
	}

	public void setZoom_factor(float zoom) {
		//check bounds
		if (zoom > 20.0f) {
			this.zoom_factor = 20.0f;
		} else if (zoom < 0.5f) {
			this.zoom_factor = 0.5f;
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
