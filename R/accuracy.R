#######################################################################
# rRDP - Interfaces RDP
# Copyright (C) 2014 Michael Hahsler and Anurag Nagar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

confusionTable <- function(actual, predicted, rank)
{
  ranka <- colnames(actual)[grep(rank, colnames(actual), ignore.case=TRUE)]
  rankp <- colnames(predicted)[grep(rank, colnames(predicted), ignore.case=TRUE)]
  actual <- factor(actual[,ranka])
  predicted <- factor(predicted[,rankp], exclude=NULL)
  
  commonLevels <- sort(unique(c(levels(actual), levels(predicted))), na.last = TRUE)
  actual <- factor(actual, levels = commonLevels, exclude=NULL)
  predicted <- factor(predicted, levels = commonLevels, exclude=NULL)
  
  table(actual,predicted, dnn=list("actual","predicted"))
}

accuracy <- function(actual, predicted, rank) {
  ct <- confusionTable(actual, predicted, rank)
  sum(diag(ct))/sum(ct)
}