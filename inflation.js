//
// Inflation demo
//
// by Andrew Pontzen 2013
//
// www.cosmocrunch.co.uk
//
// Released under Creative Commons licence BY-SA 2.5
// http://creativecommons.org/licenses/by-sa/2.5/


// will be filled with a dictionary describing the equation of state at different times:
var eos_transitions;

// how many samples to use for the plots:
var n_elements = 500;

// these arrays will have the latest plot data in them:
var a_ar = new Array(n_elements);
var adot = new Array(n_elements);
var t = new Array(n_elements);
var adot_reference = new Array(n_elements);

// the renormalization of adot, set by the boundary matching:
var adot_norm = 1.; 

// these will point to our rickshaw plots
var H_of_a, a_of_t, H_of_t, da_by_dt;


function fill_adot_end() {
    // Work backwards from a=1, where adot=1 by definition, to calculate
    // what adot was at the boundary between each "epoch".

    var adot_last = 1.0;
    var a_last = 1.0;
    var eos_last = 0.0;
    for(var i=eos_transitions.length-1; i>=0; i--) {
	var eos_i = eos_transitions[i];
	eos_i.adot_end = adot_last * Math.pow(eos_i.a_end/a_last, eos_to_Ha(eos_last)+1 );
	adot_last = eos_i.adot_end;
	eos_last = eos_i.eos;
	a_last = eos_i.a_end;
    }
}

function setup_eos_dict(a_start, a_end) {
    // Set up the different epochs

    eos_transitions = [	{ a_end: a_start, eos: 1./3, adot_end: 0 },
			{ a_end: a_end, eos: -1., adot_end: 0 },
			{ a_end: 1.0, eos: 1./3, adot_end: 0 } ]
    fill_adot_end();
}
 


function get_eos_dict(a) {
    // Get the correct desription of the equation of state for
    // the specified scalefactor.

    var a_i=0;
    var eos_i=0;
    var i=0;
    while(a_i<a) {
	a_i = eos_transitions[i].a_end;
	eos_i = eos_transitions[i];
	i=i+1;
    }
    return eos_i;
}

function eos_to_Ha(eos) {
    // Convert the equation of state to the power n of a in the
    // Friedmann integral equation, H \propto a^n
    return -(3+3*eos)/2;
}


function adot_from_a(a) {
    // Given a scalefactor a, calculate the current expansion rate
    // adot.
    var edict = get_eos_dict(a);
    return edict.adot_end*Math.pow(a/edict.a_end, eos_to_Ha(edict.eos)+1)*adot_norm;
}

function one_over_adot(a,t) {
    // Given a scalefactor a, calculate one over the current expansion rate
    // adot. Takes an additional argument t for compatibility with ODE integrator.
    return 1./adot_from_a(a);
}


function make_plot(element, x, y, log_y, ymin, ymax) {
    // Make a plot with two series, one to represent the inflationary
    // universe and one the non-inflationary universe

    var plot_data = []
    if(log_y) {
	ymin = Math.log(ymin)/Math.LN10;
	ymax = Math.log(ymax)/Math.LN10;
    }
    for(var i=1; i<n_elements; i++) {
	if(log_y) {
	    plot_data.push({x: x[i], y: Math.log(y[i])/Math.LN10});
	
	   
	} else {
	    plot_data.push({x: x[i], y: y[i]});
	}
    }
    
    graph = new Rickshaw.Graph( {
    element: document.querySelector("#"+element), 
    renderer: 'line',
    interpolation: 'linear',
    stroke: true,
    preserve: true,
	min: ymin,
	max: ymax,
    series: [{
        color: 'steelblue',
        data: plot_data,
	log_y: log_y
    }, {
        color: 'red',
        data: plot_data,
	log_y: log_y
    }]
    });

    if(log_y) {
	var yAxis = new Rickshaw.Graph.Axis.Y.Scaled( {
	    graph: graph,
	    ticks: 10,
	    scale: d3.scale.log().domain([Math.pow(10.,ymin),Math.pow(10.,ymax)])
	} ); 
	yAxis.tickValues=[1,2,5,10,20,50,100,200,500];
    } else {
	var yAxis = new Rickshaw.Graph.Axis.Y( {
	    graph: graph,
	    tickFormat: function(y) { return y.toFixed(1); }
	} ); 

    }
	
    var xAxis = new Rickshaw.Graph.Axis.X( {
        graph: graph,
	tickFormat: function(n) {
	    return n.toFixed(1);
	}
     } );
    graph.render();

    return graph;
}

function update_plot(rgraph, x, y) {
    // update an inflationary universe plot

    var plot_data = [];
    var log_y = rgraph.series[0].log_y;
    
    if(!log_y) {
	for(var i=1; i<n_elements; i++) {
	    plot_data.push({x: x[i], y: y[i]});
	}
    } else {
	for(var i=1; i<n_elements; i++) {
	    plot_data.push({x: x[i], y: Math.log(y[i])/Math.LN10});
	}
    }
    rgraph.series[0].data = plot_data;
    rgraph.update();
    rgraph.series[1].path.style["strokeDasharray"] = "6px,9px";
}



function update_visibility(a,b,c) {
    // show any warnings about the range of inflation
    if(a<=0.001 && b>=0.995) {
	$( "#warning-eternal" ).hide('fast');
	$( "#warning-dark-energy" ).hide('fast');
	$( "#warning-desitter" ).show('fast');
    } else if(a<=0.001) {
	$( "#warning-eternal" ).show('fast');
	$( "#warning-dark-energy" ).hide('fast');
	$( "#warning-desitter" ).hide('fast');	
    } else if (b>=0.995) {
	$( "#warning-eternal" ).hide('fast');
	$( "#warning-dark-energy" ).show('fast');
	$( "#warning-desitter" ).hide('fast');
    } else {
	$( "#warning-eternal" ).hide('fast');
	$( "#warning-dark-energy" ).hide('fast');
	$( "#warning-desitter" ).hide('fast');
    }	
    
    // update the commentary
    if(c>=b) {
	$("#sensible-matching").css('display','inline');
	$("#during-inflation-matching").hide();
	$(".explain_good").css('display','inline');
	$(".explain_bad").hide();
	$(".explain_intermediate").not(".explain_bad").hide();
    } else {
	$("#sensible-matching").hide();
	$("#during-inflation-matching").css('display','inline');
	if(c<a) {
	    $(".explain_bad").css('display','inline');
	    $(".explain_intermediate").not(".explain_bad").hide();
	} else {
	    $(".explain_bad").not(".explain_intermediate").hide();
	    $(".explain_intermediate").css('display','inline');
	}
	$(".explain_good").hide();
    }


    // insert the correct scalefactors into the text
    $("#a_start").html(a.toString());
    $("#a_end").html(b.toString());
    $(".a_match").html(c.toString());
    
}

function setup_a_and_adot(pivot) {
    adot_norm = 1.;

    for(var i=n_elements-1; i>=0; i--) {
	var a_this = i/n_elements;

	a_ar[i] = a_this;
	adot[i] = adot_from_a(a_ar[i]);
    }

    // now set H=1 at a=pivot
    var pivot_i = Math.round(pivot*(n_elements-1))
    if(pivot_i==0) pivot_i=1; // prevent trying to normalize on singularity

    adot_norm = adot_reference[pivot_i]/adot[pivot_i]
    if (adot_norm!=adot_norm) adot_norm=1;
    adot = numeric.mul(adot,adot_norm);
}

function make_plots() {
    // initialize our equation of state dictionary
    setup_eos_dict(0,0);

    // set up the a and adot arrays
    setup_a_and_adot(1.0);

    // calculate Hubble rate (for plotting purposes)
    var H = numeric.div(adot,a_ar);
    

    // calculate time array
    var t_sol = numeric.dopri(1.e-4,1,1.e-4,one_over_adot);
    var t_max = t_sol.at(1);
    for(var i = 0; i<n_elements; i++) {
	// store current adot for future reference
	adot_reference[i] = adot[i];
	// interpolate and store t
	t[i] = t_sol.at(a_ar[i]) - t_max;
    }


    // make the plots
    H_of_a = make_plot("H_of_a",a_ar,H,true,0.9,1e3);
    a_of_t = make_plot("a_of_t",t,a_ar,false,-0.05,1.05);
    H_of_t = make_plot("H_of_t",t,H,true,0.9,1e3);
    da_by_dt = make_plot("da_by_dt",t,adot,false,-1,11);
}

function update_plots(event,ui) {
    // update all plots for a universe with inflation spanning the
    // specified range.
    var rvals = $('#slider-range').slider('option','values');
    var inflation_start = rvals[0]/1000.;
    var inflation_end = rvals[1]/1000.;
    var matching = $('#slider-pivot').slider('option','value')/1000.
    // re-initialize our equation of state dictionary
    setup_eos_dict(inflation_start, inflation_end);

    // re-calculate adot
    setup_a_and_adot(matching);
       
    // re-calculate time array

    var t_sol = numeric.dopri(1.e-4,1,1.e-4,one_over_adot);
    var t_max = t_sol.at(1);
    for(var i = 0; i<n_elements; i++) {
	t[i] = t_sol.at(a_ar[i]) - t_max;
    }

    // kludge: if a_0=0, the numerics dictate what the range of the
    // graph is (since actually t(0)=-infinity). This can give
    // misleading effects when adjusting a_1 or a_match. So, make an
    // arbitrary t[1] choice here and fix it.
    
    if (inflation_start==0.0) {
	t[1]=-3.5;
    }

    // calculate H

    var H = numeric.div(adot,a_ar)

    // update plots

    update_plot(H_of_a,a_ar,H,true);
    update_plot(a_of_t,t,a_ar);
    update_plot(H_of_t,t,H);
    update_plot(da_by_dt,t,adot);
    update_visibility(inflation_start,inflation_end,matching);
}


$(document).ready(function () {

    // SET UP THE PAGE

    // hide the bits of text that are optional
    $('[id$="_details"]').each(function() {
	$("#"+name).hide();
	var name = this.id;
	$("#"+name+"_switch").click(function() {
	    $("#"+name).show('fast');
            $("#"+name+"_invitation").hide('fast');
	    event.preventDefault();
	    return false;
	});
    })

    $("#matching_details").hide();
    $("#timing_details").hide();
    $("#eos_details").hide();

    $("#javascript_warning").hide();
    $("#javascript_ok").show();

    // set up the jquery slider

    $( "#slider-range" ).slider({
	range: true,
	min: 0.,
	max: 1000.,
	values: [ 200., 500.], 
	slide: update_plots,
	change: update_plots,
    });

    $( "#slider-pivot" ).slider({
	min: 0.,
	max: 1000.,
	value: 1000., 
	slide: update_plots,
	change: update_plots
    });

    // set up the jquery button set graph selector

    $( "#tabs").tabs();

    // make the plots
    make_plots();

    // update the plots with a default inflationary model
    update_plots(0.2,0.5);

    // make only one plot visible
    update_visibility(0.2,0.5,1.0);

});