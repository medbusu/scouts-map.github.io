var SVY21 = function () {
  // Ref: http://www.linz.govt.nz/geodetic/conversion-coordinates/projection-conversions/transverse-mercator-preliminary-computations/index.aspx

  // WGS84 Datum
  this.a = 6378137;
  this.f = 1 / 298.257223563;

  // SVY21 Projection
  // Fundamental point: Base 7 at Pierce Resevoir.
  // Latitude: 1 22 02.9154 N, longitude: 103 49 31.9752 E (of Greenwich).

  // Known Issue: Setting (oLat, oLon) to the exact coordinates specified above
  // results in computation being slightly off. The values below give the most
  // accurate represenation of test data.
  this.oLat = 1.366666; // origin's lat in degrees
  this.oLon = 103.833333; // origin's lon in degrees
  this.oN = 38744.572; // false Northing
  this.oE = 28001.642; // false Easting
  this.k = 1; // scale factor

  this.init = function () {
    this.b = this.a * (1 - this.f);
    this.e2 = 2 * this.f - this.f * this.f;
    this.e4 = this.e2 * this.e2;
    this.e6 = this.e4 * this.e2;
    this.A0 = 1 - this.e2 / 4 - (3 * this.e4) / 64 - (5 * this.e6) / 256;
    this.A2 = (3 / 8) * (this.e2 + this.e4 / 4 + (15 * this.e6) / 128);
    this.A4 = (15 / 256) * (this.e4 + (3 * this.e6) / 4);
    this.A6 = (35 * this.e6) / 3072;
  };
  this.init();

  this.computeSVY21 = function (lat, lon) {
    //Returns a pair (N, E) representing Northings and Eastings in SVY21.

    var latR = (lat * Math.PI) / 180;
    var sinLat = Math.sin(latR);
    var sin2Lat = sinLat * sinLat;
    var cosLat = Math.cos(latR);
    var cos2Lat = cosLat * cosLat;
    var cos3Lat = cos2Lat * cosLat;
    var cos4Lat = cos3Lat * cosLat;
    var cos5Lat = cos4Lat * cosLat;
    var cos6Lat = cos5Lat * cosLat;
    var cos7Lat = cos6Lat * cosLat;

    var rho = this.calcRho(sin2Lat);
    var v = this.calcV(sin2Lat);
    var psi = v / rho;
    var t = Math.tan(latR);
    var w = ((lon - this.oLon) * Math.PI) / 180;

    var M = this.calcM(lat);
    var Mo = this.calcM(this.oLat);

    var w2 = w * w;
    var w4 = w2 * w2;
    var w6 = w4 * w2;
    var w8 = w6 * w2;

    var psi2 = psi * psi;
    var psi3 = psi2 * psi;
    var psi4 = psi3 * psi;

    var t2 = t * t;
    var t4 = t2 * t2;
    var t6 = t4 * t2;

    //	Compute Northing
    var nTerm1 = (w2 / 2) * v * sinLat * cosLat;
    var nTerm2 = (w4 / 24) * v * sinLat * cos3Lat * (4 * psi2 + psi - t2);
    var nTerm3 =
      (w6 / 720) *
      v *
      sinLat *
      cos5Lat *
      (8 * psi4 * (11 - 24 * t2) -
        28 * psi3 * (1 - 6 * t2) +
        psi2 * (1 - 32 * t2) -
        psi * 2 * t2 +
        t4);
    var nTerm4 =
      (w8 / 40320) * v * sinLat * cos7Lat * (1385 - 3111 * t2 + 543 * t4 - t6);
    var N = this.oN + this.k * (M - Mo + nTerm1 + nTerm2 + nTerm3 + nTerm4);

    //	Compute Easting
    var eTerm1 = (w2 / 6) * cos2Lat * (psi - t2);
    var eTerm2 =
      (w4 / 120) *
      cos4Lat *
      (4 * psi3 * (1 - 6 * t2) + psi2 * (1 + 8 * t2) - psi * 2 * t2 + t4);
    var eTerm3 = (w6 / 5040) * cos6Lat * (61 - 479 * t2 + 179 * t4 - t6);
    var E = this.oE + this.k * v * w * cosLat * (1 + eTerm1 + eTerm2 + eTerm3);

    return { N: N, E: E };
  };

  this.calcM = function (lat, lon) {
    var latR = (lat * Math.PI) / 180;
    return (
      this.a *
      (this.A0 * latR -
        this.A2 * Math.sin(2 * latR) +
        this.A4 * Math.sin(4 * latR) -
        this.A6 * Math.sin(6 * latR))
    );
  };

  this.calcRho = function (sin2Lat) {
    var num = this.a * (1 - this.e2);
    var denom = Math.pow(1 - this.e2 * sin2Lat, 3 / 2);
    return num / denom;
  };

  this.calcV = function (sin2Lat) {
    var poly = 1 - this.e2 * sin2Lat;
    return this.a / Math.sqrt(poly);
  };

  this.computeLatLon = function (N, E) {
    //	Returns a pair (lat, lon) representing Latitude and Longitude.

    var Nprime = N - this.oN;
    var Mo = this.calcM(this.oLat);
    var Mprime = Mo + Nprime / this.k;
    var n = (this.a - this.b) / (this.a + this.b);
    var n2 = n * n;
    var n3 = n2 * n;
    var n4 = n2 * n2;
    var G =
      this.a *
      (1 - n) *
      (1 - n2) *
      (1 + (9 * n2) / 4 + (225 * n4) / 64) *
      (Math.PI / 180);
    var sigma = (Mprime * Math.PI) / (180 * G);

    var latPrimeT1 = ((3 * n) / 2 - (27 * n3) / 32) * Math.sin(2 * sigma);
    var latPrimeT2 = ((21 * n2) / 16 - (55 * n4) / 32) * Math.sin(4 * sigma);
    var latPrimeT3 = ((151 * n3) / 96) * Math.sin(6 * sigma);
    var latPrimeT4 = ((1097 * n4) / 512) * Math.sin(8 * sigma);
    var latPrime = sigma + latPrimeT1 + latPrimeT2 + latPrimeT3 + latPrimeT4;

    var sinLatPrime = Math.sin(latPrime);
    var sin2LatPrime = sinLatPrime * sinLatPrime;

    var rhoPrime = this.calcRho(sin2LatPrime);
    var vPrime = this.calcV(sin2LatPrime);
    var psiPrime = vPrime / rhoPrime;
    var psiPrime2 = psiPrime * psiPrime;
    var psiPrime3 = psiPrime2 * psiPrime;
    var psiPrime4 = psiPrime3 * psiPrime;
    var tPrime = Math.tan(latPrime);
    var tPrime2 = tPrime * tPrime;
    var tPrime4 = tPrime2 * tPrime2;
    var tPrime6 = tPrime4 * tPrime2;
    var Eprime = E - this.oE;
    var x = Eprime / (this.k * vPrime);
    var x2 = x * x;
    var x3 = x2 * x;
    var x5 = x3 * x2;
    var x7 = x5 * x2;

    // Compute Latitude
    var latFactor = tPrime / (this.k * rhoPrime);
    var latTerm1 = latFactor * ((Eprime * x) / 2);
    var latTerm2 =
      latFactor *
      ((Eprime * x3) / 24) *
      (-4 * psiPrime2 + 9 * psiPrime * (1 - tPrime2) + 12 * tPrime2);
    var latTerm3 =
      latFactor *
      ((Eprime * x5) / 720) *
      (8 * psiPrime4 * (11 - 24 * tPrime2) -
        12 * psiPrime3 * (21 - 71 * tPrime2) +
        15 * psiPrime2 * (15 - 98 * tPrime2 + 15 * tPrime4) +
        180 * psiPrime * (5 * tPrime2 - 3 * tPrime4) +
        360 * tPrime4);
    var latTerm4 =
      latFactor *
      ((Eprime * x7) / 40320) *
      (1385 - 3633 * tPrime2 + 4095 * tPrime4 + 1575 * tPrime6);
    var lat = latPrime - latTerm1 + latTerm2 - latTerm3 + latTerm4;

    // Compute Longitude
    var secLatPrime = 1 / Math.cos(lat);
    var lonTerm1 = x * secLatPrime;
    var lonTerm2 = ((x3 * secLatPrime) / 6) * (psiPrime + 2 * tPrime2);
    var lonTerm3 =
      ((x5 * secLatPrime) / 120) *
      (-4 * psiPrime3 * (1 - 6 * tPrime2) +
        psiPrime2 * (9 - 68 * tPrime2) +
        72 * psiPrime * tPrime2 +
        24 * tPrime4);
    var lonTerm4 =
      ((x7 * secLatPrime) / 5040) *
      (61 + 662 * tPrime2 + 1320 * tPrime4 + 720 * tPrime6);
    var lon =
      (this.oLon * Math.PI) / 180 + lonTerm1 - lonTerm2 + lonTerm3 - lonTerm4;

    return { lat: lat / (Math.PI / 180), lon: lon / (Math.PI / 180) };
  };
};

// TODO: Add set styles.

mapboxgl.accessToken =
  "pk.eyJ1IjoiaXNhZHVjayIsImEiOiJjbHY0dHVydTQwY2pmMmtsb3Q4czdhaTYzIn0.1XFd8Q4sPKz9Uz1WNhwh5w";

var cv = new SVY21();
var marker = new mapboxgl.Marker();
var markerMoving = new mapboxgl.Marker(); // Initialize marker
var suggestionsContainer = document.getElementById("suggestions");
const scaleConstant = 304030427

southwest = { N: 20000, E: 2000 };
northeast = { N: 50000, E: 51000 };
southwest_wgs = cv.computeLatLon(20000, 2000);
northeast_wgs = cv.computeLatLon(50000, 51000);

console.log("Hellozzz");
// console.log(cv.computeSVY21(southwest.lat, southwest.lon));

var bounds = [
  [southwest.lon, southwest.lat], // Southwest coordinates
  [northeast.lon, northeast.lat], // Northeast coordinates
];

var mapBounds = [
  [103.6, 1.16], // Southwest coordinates
  [104.03, 1.47], // Northeast coordinates
];

var map = new mapboxgl.Map({
  container: "map",
  style: "mapbox://styles/mapbox/outdoors-v11", // You can choose different map styles
  center: [103.8198, 1.3521], // Singapore coordinates
  zoom: 12.57, // Adjust zoom level as needed
  maxBounds: mapBounds,
});

// Add scale control
map.addControl(
  new mapboxgl.ScaleControl({
    maxWidth: 80,
    unit: "metric",
  })
);

function addDataLayer() {
  map.addSource("graticule", {
    type: "geojson",
    data: graticule,
  });
  map.addLayer({
    id: "graticule",
    type: "line",
    source: "graticule",
  });
  map.addLayer({
    id: "countour-labels",
    type: "symbol",
    source: {
      type: 'vector',
      url: 'mapbox://mapbox.mapbox-terrain-v2'
    },
    "source-layer": "contour",
    'layout': {
      'visibility': 'visible',
      'symbol-placement': 'line',
      'text-field': ['concat', ['to-string', ['get', 'ele']], 'm']
    },
    'paint': {
      'icon-color': '#877b59',
      'icon-halo-width': 1,
      'text-color': '#877b59',
      'text-halo-width': 1
    },
    'minzoom': 4
  });

  
  
  // map.addLayer({
  //   "id": "countours",
  //   "type": "line",
  //   "source": {
  //     type: 'vector',
  //     url: 'mapbox://mapbox.mapbox-terrain-v2'
  //   },
  //   "source-layer": "contour",
  //   'layout': {
  //     'visibility': 'visible',
  //     'line-join': 'round',
  //     'line-cap': 'round'
  //   },
  //   'paint': {
  //     'line-color': '#877b59',
  //     'line-width': 1
  //   }
  // })
}
map.on("styledata", function () {
  // Triggered when `setStyle` is called.
  addDataLayer();
});

map.on("load", function () {
  // Generate grid squares
  // generateGridSquares();
  addDataLayer();

  
});

var crosshair = document.getElementById("crosshair");

map.on("mousemove", function (e) {
  // Update the position of the crosshair based on the mouse position
  crosshair.style.left = e.point.x + "px";
  crosshair.style.top = e.point.y + "px";
});

const graticule = {
  type: "FeatureCollection",
  features: [],
};

for (let easting = southwest.E; easting <= northeast.E; easting += 1000) {
  let svy21 = cv.computeLatLon(southwest.N, easting);
  // console.log("southwest is", southwest.lat, " + ", northeast.lat);

  graticule.features.push({
    type: "Feature",
    geometry: {
      type: "LineString",
      coordinates: [
        [svy21.lon, southwest_wgs.lat],
        [svy21.lon, northeast_wgs.lat],
      ],
    },
    properties: { value: easting },
  });
}
for (let northing = southwest.N; northing <= northeast.N; northing += 1000) {
  let svy21 = cv.computeLatLon(northing, southwest.E);
  graticule.features.push({
    type: "Feature",
    geometry: {
      type: "LineString",
      coordinates: [
        [southwest_wgs.lon, svy21.lat],
        [northeast_wgs.lon, svy21.lat],
      ],
    },
    properties: { value: northing },
  });
}

map.on("click", function (e) {
  var svy21Coords = cv.computeSVY21(e.lngLat.lat, e.lngLat.lng);
  document.getElementById("coordinates").innerText =
    "SVY21 Coordinates: " +
    svy21Coords.N.toFixed(2) +
    ", " +
    svy21Coords.E.toFixed(2);

  marker.setLngLat(e.lngLat).addTo(map);
});

document.getElementById("search").addEventListener("input", function (e) {
  var query = e.target.value;
  console.log(query);
  if (query.trim() !== "") {
    // Perform geocoding API request
    fetch(
      "https://api.mapbox.com/geocoding/v5/mapbox.places/" +
        encodeURIComponent(query) +
        ".json?access_token=" +
        mapboxgl.accessToken
    )
      .then((response) => response.json())
      .then((data) => {
        showSuggestions(data.features);
      })
      .catch((error) => {
        console.error("Error:", error);
      });
    console.log("hello");
  } else {
    console.log("no input");
    clearSuggestions();
  }
});

function showSuggestions(features) {
  clearSuggestions();
  console.log(features);
  features.forEach((feature) => {
    var suggestion = document.createElement("div");
    console.log(suggestion);
    suggestion.classList.add("suggestion");
    suggestion.textContent = feature.place_name;
    suggestion.addEventListener("click", function () {
      var coordinates = feature.center;
      // Move the map to the searched location
      map.flyTo({ center: coordinates, zoom: 12.57 });
      var svy21Coords = cv.computeSVY21(coordinates[1], coordinates[0]);
      console.log(svy21Coords);
      // Update textbox value with SVY21 coordinates
      document.getElementById("coordinates").innerText =
        "SVY21 Coordinates: " +
        svy21Coords.N.toFixed(2) +
        ", " +
        svy21Coords.E.toFixed(2);

      // Move the marker to the searched location
      marker.setLngLat(coordinates).addTo(map);
      // Clear suggestions
      clearSuggestions();
      // Clear search input
      document.getElementById("search").value = "";
    });
    suggestionsContainer.appendChild(suggestion);
  });
}

function clearSuggestions() {
  suggestionsContainer.innerHTML = "";
}

// Add event listener to the button
const layerList = document.getElementById("menu");
const inputs = layerList.getElementsByTagName("input");

for (const input of inputs) {
  input.onclick = (layer) => {
    const layerId = layer.target.id;
    map.setStyle("mapbox://styles/mapbox/" + layerId);
  };
}

document.getElementById("setScaleButton").addEventListener("click", function () {
    // Calculate the zoom level for a scale of 1:50000
    var zoomLevel = 12.57
    console.log("Hi", zoomLevel)

    // Set the map's zoom level
    map.setZoom(zoomLevel);
  });

map.on('zoom', function () {
  var zoomLevel = document.getElementById('zoom-level');
  var scaleLevel = document.getElementById('scale-level');
  var level = map.getZoom();
  var scaleLevelDisplay = scaleConstant / (2**level);
  zoomLevel.innerText = 'Scale Bar = 1:' + scaleLevelDisplay.toFixed(0);
  scaleLevel.innerText = 'Zoom Level: ' + level.toFixed(3);
  // console.log("LOL", zoomLevel.innerText);
  // 'Zoom Level: ' + level.toFixed(3) + 
});

document.getElementById("in").addEventListener("click", function () {
  // Calculate the zoom level for a scale of 1:50000
  var zoomLevel = map.getZoom()
  // Set the map's zoom level
  map.setZoom(zoomLevel + 1);
});

document.getElementById("out").addEventListener("click", function () {
  // Calculate the zoom level for a scale of 1:50000
  var zoomLevel = map.getZoom()
  // Set the map's zoom level
  map.setZoom(zoomLevel - 1);
});



function setZoomLevel() {
  var zoomLevel = parseFloat(document.getElementById('zoomLevelInput').value);
  if (isNaN(zoomLevel) || zoomLevel < 0) {
      alert('Please enter a valid positive number for the zoom level.');
      return;
  }
  
  let zoomLevelAdj = Math.log2(scaleConstant/zoomLevel)
  map.setZoom(zoomLevelAdj);
}