/**
* @file flan.h
* @brief Header file for flan.cpp.
*/
#pragma once

#include "flan_types.h" // For Inputs type

/**
* @brief Entry point to Flan library
*
* @param inpts User defined input options that, along with the default options,
*        control a Flan simulation.
*/
void flan(const Inputs& inpts);
