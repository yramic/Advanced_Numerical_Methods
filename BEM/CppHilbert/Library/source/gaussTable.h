///////////////////////////////////////////////////////////////////////////////
/// \file gaussTable.h
/// \brief This file is part of the HILBERT program package for the numerical
///        solution of the Laplace equation with mixed boundary conditions by
///        use of BEM in 2D.
///
///  It provides gauss points and weights for order 2,4,8,16, and 32.
///  Furthermore, points and weights for a 7-point Gaussian quadrature rule on
///  triangles is provided. The data should not be accessed directly. Use the
///  functions declared in gaussQuadrature.h instead.
///
/// VERSION: 3.1
/// (C) 2009-2013 HILBERT-Team '09, '10, '12
/// support + bug report:  hilbert@asc.tuwien.ac.at
///////////////////////////////////////////////////////////////////////////////
#ifndef _GAUSS_TABLE_H_GUARD_
#define _GAUSS_TABLE_H_GUARD_

/* 7-point quadrature on Reference Triangle conv{(0,0),(1,0),(0,1)} */
const static double gaussT7_x[] = {
    0.101286507323456, 0.797426985353087, 0.101286507323456, 0.470142064105115,
    0.470142064105115, 0.059715871789770, 0.333333333333333};

const static double gaussT7_y[] = {
    0.101286507323456, 0.101286507323456, 0.797426985353087, 0.059715871789770,
    0.470142064105115, 0.470142064105115, 0.333333333333333};

const static double gaussT7_w[] = {0.062969590272414,
                                   0.062969590272414,
                                   0.062969590272414,
                                   0.066197076394253,
                                   0.066197076394253,
                                   0.066197076394253,
                                   0.1125};

/* 2-point Gauss quadrature: */
const static double gauss2_points[] = {-0.5773502691896257645091488,
                                       0.5773502691896257645091488};

const static double gauss2_weights[] = {1, 1};

/* 4-point Gauss quadrature:
 */
const static double gauss4_points[] = {
    -0.8611363115940525752239465, -0.3399810435848562648026658,
    0.3399810435848562648026658, 0.8611363115940525752239465};

const static double gauss4_weights[] = {
    0.3478548451374538573730639, 0.6521451548625461426269361,
    0.6521451548625461426269361, 0.3478548451374538573730639};

/* 8-point Gauss quadrature:
 */
const static double gauss8_points[] = {
    -0.9602898564975362316835609, -0.7966664774136267395915539,
    -0.5255324099163289858177390, -0.1834346424956498049394761,
    0.1834346424956498049394761,  0.5255324099163289858177390,
    0.7966664774136267395915539,  0.9602898564975362316835609};

const static double gauss8_weights[] = {
    0.1012285362903762591525314, 0.2223810344533744705443560,
    0.3137066458778872873379622, 0.3626837833783619829651504,
    0.3626837833783619829651504, 0.3137066458778872873379622,
    0.2223810344533744705443560, 0.1012285362903762591525314};

/* 16-point Gauss quadrature:
 */
const static double gauss16_points[] = {
    -0.9894009349916499325961542, -0.9445750230732325760779884,
    -0.8656312023878317438804679, -0.7554044083550030338951012,
    -0.6178762444026437484466718, -0.4580167776572273863424194,
    -0.2816035507792589132304605, -0.0950125098376374401853193,
    0.0950125098376374401853193,  0.2816035507792589132304605,
    0.4580167776572273863424194,  0.6178762444026437484466718,
    0.7554044083550030338951012,  0.8656312023878317438804679,
    0.9445750230732325760779884,  0.9894009349916499325961542};

const static double gauss16_weights[] = {
    0.0271524594117540948517806, 0.0622535239386478928628438,
    0.0951585116824927848099251, 0.1246289712555338720524763,
    0.1495959888165767320815017, 0.1691565193950025381893121,
    0.1826034150449235888667637, 0.1894506104550684962853967,
    0.1894506104550684962853967, 0.1826034150449235888667637,
    0.1691565193950025381893121, 0.1495959888165767320815017,
    0.1246289712555338720524763, 0.0951585116824927848099251,
    0.0622535239386478928628438, 0.0271524594117540948517806};

/* 32-point Gauss quadrature
 */
const static double gauss32_points[] = {
    -0.9972638618494815635449811, -0.9856115115452683354001750,
    -0.9647622555875064307738119, -0.9349060759377396891709191,
    -0.8963211557660521239653072, -0.8493676137325699701336930,
    -0.7944837959679424069630973, -0.7321821187402896803874267,
    -0.6630442669302152009751152, -0.5877157572407623290407455,
    -0.5068999089322293900237475, -0.4213512761306353453641194,
    -0.3318686022821276497799168, -0.2392873622521370745446032,
    -0.1444719615827964934851864, -0.0483076656877383162348126,
    0.0483076656877383162348126,  0.1444719615827964934851864,
    0.2392873622521370745446032,  0.3318686022821276497799168,
    0.4213512761306353453641194,  0.5068999089322293900237475,
    0.5877157572407623290407455,  0.6630442669302152009751152,
    0.7321821187402896803874267,  0.7944837959679424069630973,
    0.8493676137325699701336930,  0.8963211557660521239653072,
    0.9349060759377396891709191,  0.9647622555875064307738119,
    0.9856115115452683354001750,  0.9972638618494815635449811};

const static double gauss32_weights[] = {
    0.0070186100094700966004071, 0.0162743947309056706051706,
    0.0253920653092620594557526, 0.0342738629130214331026877,
    0.0428358980222266806568786, 0.0509980592623761761961632,
    0.0586840934785355471452836, 0.0658222227763618468376501,
    0.0723457941088485062253994, 0.0781938957870703064717409,
    0.0833119242269467552221991, 0.0876520930044038111427715,
    0.0911738786957638847128686, 0.0938443990808045656391802,
    0.0956387200792748594190820, 0.0965400885147278005667648,
    0.0965400885147278005667648, 0.0956387200792748594190820,
    0.0938443990808045656391802, 0.0911738786957638847128686,
    0.0876520930044038111427715, 0.0833119242269467552221991,
    0.0781938957870703064717409, 0.0723457941088485062253994,
    0.0658222227763618468376501, 0.0586840934785355471452836,
    0.0509980592623761761961632, 0.0428358980222266806568786,
    0.0342738629130214331026877, 0.0253920653092620594557526,
    0.0162743947309056706051706, 0.0070186100094700966004071};

#endif
