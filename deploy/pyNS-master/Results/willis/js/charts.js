function render_charts(mesh) {
	if($.browser.msie){
	ie_rendering_charts(mesh);
	}
	else {
	rendering_charts(mesh);
	}
}

function render_charts_adapt(mesh) {
	if($.browser.msie){
	ie_rendering_charts_adapt(mesh);
	}
	else {
	rendering_charts_adapt(mesh);
	}
}

function ie_rendering_charts(mesh) {
    $('#mesh_info').removeClass("hide");
    jQuery.getJSON('json/'+mesh+'.json', function(data) {
    	$('[id^="mesh_"]').removeClass("active");
    	$('[id^="mesh_"]').children('a').removeClass("selected");
    	$('#mesh_'+data.name).addClass("active");
    	$('#mesh_'+data.name).children('a').addClass("selected");
    	$('#length').html(data.length);
     	$('#min_d').html(data.diameter_min);
     	$('#max_d').html(data.diameter_max);
     	$('#mean_p').html(data.mean_pressure);
     	$('#mean_q').html(data.mean_flow);
     	$('#mean_w').html(data.mean_wss);
     	$('#mean_r').html(data.mean_re);
     	
    	$.each(data.items, function(i,item){
    	
    	var chart1 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "flow",
                    type: 'line',
					height: 350
                },
    			title: {
                    text: 'Volumetric Flow Rate'
                },
    			xAxis: {
                    title: {
                        text: 'Time (s)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Volumetric Flow Rate (mL/min)'
                    },
                    min: data.min_flow
                },
                plotOptions: {
                	series: {
                	lineWidth:4,
                	marker: {
                    	enabled: false,
                    	states: {
                        	hover: {
                            	enabled: true
                        	}
                    	}
               		  }
                	}
                },
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Flow volume", data:item.flow}]  
            }  
        );    
        
        var chart2 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "pressure",
                    type: 'line',
					height: 350
                },
    			title: {
                    text: 'Pressure'
                },
    			xAxis: {
                    title: {
                        text: 'Time (s)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Pressure signal (mmHg)'
                    },
                    min : data.min_pressure
                },
                plotOptions: {
                	series: {
                	lineWidth:4,
                	marker: {
                    	enabled: false,
                    	states: {
                        	hover: {
                            	enabled: true
                        	}
                    	}
               		  }
                	}
                },
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Pressure", data:item.pressure}]  
            }  
        ); 
        
        var chart3 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "wss",
                    type: 'line',
					height: 350
                },
    			title: {
                    text: 'Wall Shear Stress'
                },
    			xAxis: {
                    title: {
                        text: 'Time (s)'
                    }
                },
    			yAxis: {
                    title: {
                		text: 'Wall Shear Stress (dynes/cm^2)',
                		enabled:true
           			},
                    min : data.min_wss
                },
                plotOptions: {
                	series: {
                	lineWidth:4,
                	marker: {
                    	enabled: false,
                    	states: {
                        	hover: {
                            	enabled: true
                        	}
                    	}
               		  }
                	}
                },
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Wss", data:item.wss}]  
            }  
        ); 
        
        var chart4 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "re",
                    type: 'line',
					height: 350
                },
    			title: {
                    text: 'Reynolds Number'
                    },
    			xAxis: {
                    title: {
                        text: 'Time (s)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Reynolds Number'
                    },
                    min : data.min_re
                },
                plotOptions: {
                	series: {
                	lineWidth:4,
                	marker: {
                    	enabled: false,
                    	states: {
                        	hover: {
                            	enabled: true
                        	}
                    	}
               		  }
                	}
                },
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Re", data:item.re}]  
            }  
        ); 
        
        });
    });
    }
	
function rendering_charts(mesh) { 
    $('#mesh_info').removeClass("hide");
   
    jQuery.getJSON('json/'+mesh+'.json', function(data) {
    	
    	$('[id^="mesh_"]').removeClass("active");
    	$('[id^="mesh_"]').children('a').removeClass("selected");
    	$('#mesh_'+data.name).addClass("active");
    	$('#mesh_'+data.name).children('a').addClass("selected");
    	
    	$('#length').html(data.length);
     	$('#min_d').html(data.diameter_min);
     	$('#max_d').html(data.diameter_max);
     	$('#mean_p').html(data.mean_pressure);
     	$('#mean_q').html(data.mean_flow);
     	$('#mean_w').html(data.mean_wss);
     	$('#mean_r').html(data.mean_re);
     	
    	$.each(data.items, function(i,item){
    	
    	var chart1 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "flow",
                    type: 'line'
                },
    			title: {
                    text: 'Volumetric Flow Rate'
                },
    			xAxis: {
                    title: {
                        text: 'Time (s)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Volumetric Flow Rate (mL/min)'
                    },
                    min: data.min_flow
                },
                plotOptions: {
                	series: {
                	lineWidth:4,
                	marker: {
                    	enabled: false,
                    	states: {
                        	hover: {
                            	enabled: true
                        	}
                    	}
               		  }
                	}
                },
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Flow volume", data:item.flow}]  
            }  
        );    
        
        var chart2 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "pressure",
                    type: 'line'
                },
    			title: {
                    text: 'Pressure'
                },
    			xAxis: {
                    title: {
                        text: 'Time (s)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Pressure signal (mmHg)'
                    },
                    min : data.min_pressure
                },
                plotOptions: {
                	series: {
                	lineWidth:4,
                	marker: {
                    	enabled: false,
                    	states: {
                        	hover: {
                            	enabled: true
                        	}
                    	}
               		  }
                	}
                },
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Pressure", data:item.pressure}]  
            }  
        ); 
        
        var chart3 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "wss",
                    type: 'line'
                },
    			title: {
                    text: 'Wall Shear Stress'
                },
    			xAxis: {
                    title: {
                        text: 'Time (s)'
                    }
                },
    			yAxis: {
                    title: {
                		text: 'Wall Shear Stress (dynes/cm^2)',
                		enabled:true
           			},
                    min : data.min_wss
                },
                plotOptions: {
                	series: {
                	lineWidth:4,
                	marker: {
                    	enabled: false,
                    	states: {
                        	hover: {
                            	enabled: true
                        	}
                    	}
               		  }
                	}
                },
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Wss", data:item.wss}]  
            }  
        ); 
        
        var chart4 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "re",
                    type: 'line'
                },
    			title: {
                    text: 'Reynolds Number'
                    },
    			xAxis: {
                    title: {
                        text: 'Time (s)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Reynolds Number'
                    },
                    min : data.min_re
                },
                plotOptions: {
                	series: {
                	lineWidth:4,
                	marker: {
                    	enabled: false,
                    	states: {
                        	hover: {
                            	enabled: true
                        	}
                    	}
               		  }
                	}
                },
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Re", data:item.re}]  
            }  
        ); 
        
        });
        
    });
    
    }
    function ie_rendering_charts_adapt(mesh) { 
    jQuery.getJSON('json/'+mesh+'.json', function(data) {
    	$('[id^="mesh_"]').removeClass("active");
    	$('[id^="mesh_"]').children('a').removeClass("selected");
    	$('#mesh_'+data.name).addClass("active");
    	$('#mesh_'+data.name).children('a').addClass("selected");
		$.each(data.items, function(i,item){
    	var chart1 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "flow",
                    type: 'spline',
					height: 350
                },
    			title: {
                    text: 'Volumetric flow Rate'
                },
    			xAxis: {
                    title: {
                        text: 'Time (days)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Volumetric Flow Rate (mL/min)'
                    },
					min: data.min_q
                },
                
                plotOptions: {
					spline: {
						marker: {
							radius: 4,
							lineColor: '#666666',
							lineWidth: 1
						}
					}
				},
				
				tooltip: {
					crosshairs: true,
					shared: true
				},
                        
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Flow volume", data:item.flow}]  
            }  
        );    
        
        var chart2 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "wss",
                    type: 'spline',
					height: 350
                },
    			title: {
                    text: 'Pressure'
                },
    			xAxis: {
                    title: {
                        text: 'Time (days)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Pressure signal (mmHg)'
                    },
                    min : 0
                },
                
                plotOptions: {
					spline: {
						marker: {
							radius: 4,
							lineColor: '#666666',
							lineWidth: 1
						}
					}
				},
				
				tooltip: {
					crosshairs: true,
					shared: true
				},
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Pressure", data:item.pressure}]  
            }  
        ); 
        
        var chart3 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "re",
                    type: 'spline',
					height: 350
                },
    			title: {
                    text: 'Wall Shear Stress Peak'
                },
    			xAxis: {
                    title: {
                        text: 'Time (days)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Wall Shear Stress Peak (dynes/cm2)'
                    },
                    min : 0
                },
                
                plotOptions: {
					spline: {
						marker: {
							radius: 4,
							lineColor: '#666666',
							lineWidth: 1
						}
					}
				},
				
                tooltip: {
					crosshairs: true,
					shared: true
				},
				
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Wss Peak", data:item.wssP}]  
            }  
        ); 
        
        var chart4 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "pressure",
                    type: 'spline',
					height: 350
                },
    			title: {
                    text: 'Diameter'
                    },
    			xAxis: {
                    title: {
                        text: 'Time (days)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Diameter (mm)'
                    },
                    min : 0
                },
                
                plotOptions: {
					spline: {
						marker: {
							radius: 4,
							lineColor: '#666666',
							lineWidth: 1
						}
					}
				},
				
				tooltip: {
					crosshairs: true,
					shared: true
				},
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Diameter", data:item.diameter}]  
            }  
        ); 
        
        });
    });
    }function rendering_charts_adapt(mesh) { 
    jQuery.getJSON('json/'+mesh+'.json', function(data) {
    	$('[id^="mesh_"]').removeClass("active");
    	$('[id^="mesh_"]').children('a').removeClass("selected");
    	$('#mesh_'+data.name).addClass("active");
    	$('#mesh_'+data.name).children('a').addClass("selected");
		$.each(data.items, function(i,item){
    	var chart1 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "flow",
                    type: 'spline'
                },
    			title: {
                    text: 'Volumetric flow Rate'
                },
    			xAxis: {
                    title: {
                        text: 'Time (days)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Volumetric Flow Rate (mL/min)'
                    },
					min: data.min_q
                },
                
                plotOptions: {
					spline: {
						marker: {
							radius: 4,
							lineColor: '#666666',
							lineWidth: 1
						}
					}
				},
				
				tooltip: {
					crosshairs: true,
					shared: true
				},
                        
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Flow volume", data:item.flow}]  
            }  
        );    
        
        var chart2 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "wss",
                    type: 'spline'
                },
    			title: {
                    text: 'Pressure'
                },
    			xAxis: {
                    title: {
                        text: 'Time (days)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Pressure signal (mmHg)'
                    },
                    min : 0
                },
                
                plotOptions: {
					spline: {
						marker: {
							radius: 4,
							lineColor: '#666666',
							lineWidth: 1
						}
					}
				},
				
				tooltip: {
					crosshairs: true,
					shared: true
				},
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Pressure", data:item.pressure}]  
            }  
        ); 
        
        var chart3 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "re",
                    type: 'spline'
                },
    			title: {
                    text: 'Wall Shear Stress Peak'
                },
    			xAxis: {
                    title: {
                        text: 'Time (days)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Wall Shear Stress Peak (dynes/cm2)'
                    },
                    min : 0
                },
                
                plotOptions: {
					spline: {
						marker: {
							radius: 4,
							lineColor: '#666666',
							lineWidth: 1
						}
					}
				},
				
                tooltip: {
					crosshairs: true,
					shared: true
				},
				
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Wss Peak", data:item.wssP}]  
            }  
        ); 
        
        var chart4 = new Highcharts.Chart(
    	{
                chart:{
                    renderTo: "pressure",
                    type: 'spline'
                },
    			title: {
                    text: 'Diameter'
                    },
    			xAxis: {
                    title: {
                        text: 'Time (days)'
                    }
                },
    			yAxis: {
                    title: {
                        text: 'Diameter (mm)'
                    },
                    min : 0
                },
                
                plotOptions: {
					spline: {
						marker: {
							radius: 4,
							lineColor: '#666666',
							lineWidth: 1
						}
					}
				},
				
				tooltip: {
					crosshairs: true,
					shared: true
				},
                
                credits: {
                	enabled: false
                },
                
                legend: {
                	enabled: false
                },
                
                exporting: {
                	filename : data.name
                },
                
                series: [{name:"Diameter", data:item.diameter}]  
            }  
        ); 
        
        });
    });
    }