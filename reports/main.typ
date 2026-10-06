/////////////
// Imports //
/////////////
#import "@preview/adaptable-pset:0.2.0": *
#import "@preview/physica:0.9.8": *
#import "@preview/unify:0.8.1": *
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#show: codly-init.with()
#codly(languages: codly-languages)

/////////////////
// Maths Setup //
/////////////////

// upright vectors
#let vectorboldupright(a) = vb($upright(#a)$)
#let vbu = vectorboldupright
#let vectorunitupright(a) = vu($upright(#a)$)
#let vuu = vectorunitupright
#let vectorarrowupright(a) = va($upright(#a)$)
#let vau = vectorarrowupright

// automatically use square brackets for vectors and matricies
#set math.vec(delim: "[")
#set math.mat(delim: "[")
#let vecrowOld = vecrow
#let vecrow = vecrowOld.with(delim: "[")

////////////////////
// Document Setup //
////////////////////

// assignment info
#show: homework.with(
    title: "HW03",
    author: "Vai Srivastava",
    collaborators: [],
    course-id: "ENAE 601: Astrodynamics",
    instructor: "Dr. Healy",
    semester: "Fall 2026",
    due-time: datetime(
        // due date
        year: 2026,
        month: 10,
        day: 07,

        // due time
        hour: 23,
        minute: 59,
        second: 00,
    ).display(
        "[month repr:long] [day padding:zero], [year] at [hour repr:24]:[minute padding:zero]:[second padding:zero]",
    ),

    // (defaults to A4)
    paper-size: "us-letter",
)

// document settings
#set text(font: "New Computer Modern", size: 10pt)
#set enum(numbering: "a)")

// problem headings
#let probOld = prob
#let prob = prob.with(color: black)

////////////////////////////
// The Assignment Itself: //
// Problems and Solutions //
////////////////////////////

// NOTE: Where classical orbital elements are requested, use standard orbital elements with semimajor axis $a$---not angular momentum magnitude $h$---as per #link("HW02")[../../hw02/reports/main.typ]

#prob(title: "Curtis 4.3")[
    For a geocentric satellite whose inertial position and velocity vectors in a geocentric equatorial frame are:
    $
        vbu(r) = & 2500 vuu(I) & +16000 vuu(J) & +4000 vuu(K) &   space (unit("km")) \
        vbu(v) = & -3 vuu(I)   &       -vuu(J) & +5 vuu(K)    & space (unit("km/s"))
    $

    + Find the classical orbital elements $e$, $a$, $i$, $Omega$, $omega$, $theta$.
    <hwk:p01a>

    + Compute the equinoctial orbital elements $a$, $h$ #emph[(this is not angular momentum)], $k$, $p$ #emph[(this is not semilatus rectum)], $q$, $lambda$.
    <hwk:p01b>
] <hwk:p01>

// TODO: answer
+ Answer
<hwk:s01a>

// TODO: answer
+ Answer
<hwk:s01b>

#pagebreak(weak: true)

#prob(title: "Curtis 4.4")[
    At a given instant, the position $vbu(r)$ and velocity $vbu(v)$ of a satellite in the geocentric equatorial frame are:
    $
        vbu(r) & =          &           & -13000 vuu(K) &   space (unit("km")) \
        vbu(v) & = 4 vuu(I) & +3 vuu(J) & -3 vuu(K)     & space (unit("km/s"))
    $

    + Find the classical orbital elements $e$, $a$, $i$, $Omega$, $omega$, $theta$.
    <hwk:p02a>

    + Compute the equinoctial orbital elements $a$, $h$ #emph[(this is not angular momentum)], $k$, $p$ #emph[(this is not semilatus rectum)], $q$, $lambda$.
    <hwk:p02b>
] <hwk:p02>

// TODO: answer
+ Answer
<hwk:s02a>

// TODO: answer
+ Answer
<hwk:s02b>

#pagebreak(weak: true)

#prob(title: "Curtis 4.5")[
    At time $t_0$ (relative to perigee passage) the position $vbu(r)$ and velocity $vbu(v)$ of a satellite in the geocentric equatorial frame are:
    $
        vbu(r) & = 6500 vuu(I) & -7500 vuu(J) & -2500 vuu(K) &   space (unit("km")) \
        vbu(v) & = 4 vuu(I)    &    +3 vuu(J) & -3 vuu(K)    & space (unit("km/s"))
    $

    + Find the classical orbital elements $e$, $a$, $i$, $Omega$, $omega$, $theta$.
    <hwk:p03a>

    + Compute the equinoctial orbital elements $a$, $h$ #emph[(this is not angular momentum)], $k$, $p$ #emph[(this is not semilatus rectum)], $q$, $lambda$.
    <hwk:p03b>
] <hwk:p03>

// TODO: answer
+ Answer
<hwk:s03a>

// TODO: answer
+ Answer
<hwk:s03b>

#pagebreak(weak: true)

#prob(title: "Curtis 4.6")[
    With respect to the geocentric equatorial frame, the position vector of a spacecraft is $vbu(r) = - 6000 vuu(I) - 1000 vuu(J) - 5000 vuu(K) space (unit("km"))$ and the orbit's eccentricity vector is $vbu(e) = 0.4 vuu(I) + 0.5 vuu(J) + 0.6 vuu(K)$. Calculate the true anomaly $theta$ if the satellite is approaching perigee.
] <hwk:p04>

// TODO: answer
+ Answer
<hwk:s04>

#pagebreak(weak: true)

#prob(title: "Curtis 4.7")[
    Given that, relative to the geocentric equatorial frame, $vbu(r) = - 6600 vuu(I) - 1300 vuu(J) - 5200 vuu(K) space (unit("km"))$, the eccentricity vector is $vbu(e) = - 0.4 vuu(I) - 0.5 vuu(J) - 0.6 vuu(K)$, and the satellite is flying toward perigee, calculate the inclination of the orbit.
] <hwk:p05>

// TODO: answer
+ Answer
<hwk:s05>

#pagebreak(weak: true)

#prob(title: "")[
    Find the equinoctial elements for the following orbits:

    + Position and Velocity:
        $
            vbu(r) = mat(3698.3455; -34232.4473; 0) space (unit("km")), quad
            vbu(v) = mat(2.1482; -1.4944; 0) space (unit("km/s"))
        $
    <hwk:p06a>

    + Position and Velocity:
        $
            vbu(r) = mat(-21082.0848; 36515.242; 0) space (unit("km")), quad
            vbu(v) = mat(-2.6627; -1.5373; 0) space (unit("km/s"))
        $
    <hwk:p06b>

    + Position and Velocity:
        $
            vbu(r) = mat(-15843.4562; -2247.7766; 21200.4627) space (unit("km")), quad
            vbu(v) = mat(1.4544; -3.5188; 0.7138) space (unit("km/s"))
        $
    <hwk:p06c>
] <hwk:p06>

// TODO: answer
+ Answer
<hwk:s06a>

// TODO: answer
+ Answer
<hwk:s06b>

// TODO: answer
+ Answer
<hwk:s06c>

#pagebreak(weak: true)

== Code

#codly(header: [./src/index.py])
#raw(read("../src/index.py"), block: true, lang: "python") <code:index.py>
