// ### Extract NDVI for IRG
// # Alec Robitaille
// # November 2018

// ### Image Import
// NDVI
var ndvi = ee.ImageCollection('MODIS/006/MOD13Q1');


// ### Feature Collection Import
// Points
var reg = 'reg-pts-';
var numb = 1;
var root = 'users/alr201/48hr/';
var points = ee.FeatureCollection(root + reg + numb);

// ### Filter
ndvi = ndvi.filter(ee.Filter.calendarRange(2007, 2013, 'year'));

// ### Add year
var addYear = function(img) {
  return(img.addBands(ee.Image(img.date().get('year')).rename('yr')));
};

ndvi = ndvi.map(addYear);

// ### Sample
var samplePts = function(img) {
  return img.reduceRegions(points, ee.Reducer.mean())
            .copyProperties(points);
};

var reduced = ndvi.map(samplePts)
									.flatten()
									.filter(ee.Filter.neq('NDVI', null));
// ### Export
var outputString  = reg + numb;
Export.table.toDrive(reduced, outputString, 'IRG');
