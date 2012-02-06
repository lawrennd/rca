% function includeMotionCategories()

%{
% % Categories, [1-9]
% HUMAN_INTERACTION = 1;
% INTERACTION_WITH_ENVIRONMENT = 2;
% LOCOMOTION = 3;
% PHYSICAL_ACTIVITIES_AND_SPORTS = 4;
% SITUATIONS_AND_SCENARIOS = 5;
% 
% % Sub-categories [10-99]
% TWO_SUBJECTS = 10;
% PLAYGROUND = 20;
% UNEVEN_TERRAIN = 21;
% PATH_WITH_OBSTACLES = 22;
% RUNNING = 30;
% WALKING = 31;
% JUMPING = 32;
% VARIED = 33;
% BASKETBALL = 40;
% DANCE = 41;
% GYMNASTICS = 42;
% ACROBATICS = 43;
% MARTIAL_ARTS = 44;
% SOCCER = 45;
% BOXING = 46;
% GENERAL_EXERCISE_AND_STRETCHING = 47;
% COMMON_BEHAVIORS_AND_EXPRESSIONS = 50;
% PANTOMIME = 51;
% COMMUNICATION_GESTURES_AND_SIGNALS = 52;
% CROSS_CATEGORY_VARIETY = 53;
%
% a = struct('CategoryName', HUMAN_INTERACTION, 'Subcategories', TWO_SUBJECTS);
% b = struct('CategoryName', INTERACTION_WITH_ENVIRONMENT, 'Subcategories', [PLAYGROUND, UNEVEN_TERRAIN, PATH_WITH_OBSTACLES]);
% c = struct('CategoryName', LOCOMOTION, 'Subcategories', [RUNNING, WALKING, JUMPING, VARIED]);
% d = struct('CategoryName', PHYSICAL_ACTIVITIES_AND_SPORTS, 'Subcategories', ...
%     [BASKETBALL, DANCE, GYMNASTICS, ACROBATICS, MARTIAL_ARTS, SOCCER, BOXING, GENERAL_EXERCISE_AND_STRETCHING]);
% e = struct('CategoryName', SITUATIONS_AND_SCENARIOS, 'Subcategories', [COMMON_BEHAVIORS_AND_EXPRESSIONS, PANTOMIME, COMMUNICATION_GESTURES_AND_SIGNALS, CROSS_CATEGORY_VARIETY]);
%}

% Sub-category structs
two_subjects = struct( 'name', 'two_subjects', 'subject', {18,19,20,21,22,23,33,34} , 'trial', {1:15, 1:15, 1:13, 1:13, 1:25, 1:25, 1:2, 1:2});
playground = struct( 'name', 'playground', 'subject', {1,43} , 'trial', {1:14, 2:3});
uneven_terrain = struct( 'name', 'uneven_terrain' , 'subject', {3,36} , 'trial', {1:4, [1,4:8,10:37]});
path_with_obstacles = struct( 'name', 'path_with_obstacles', 'subject', {35,54} , 'trial', {27, 20:22});
running = struct( 'name', 'running' , 'subject', {2,9,16,35,38} , 'trial', {3, 1:11, [8,35:46,48:57], 17:26, 3});
walking = struct( 'name', 'walking', 'subject', {2,5,6,7,8,9,10,12,15,16,17,26,27,29,32,35,36,37,38,39,40,41,43,45,46,47,49,55,56} , ...
    'trial', {1:2, 1, 1, 1:12, 1:11, 12, 4, 1:3, [1,3,9,14], [11:34,47,58], 1:9, 1, 1, 1, 1:2, [1:16,28:34], [2,3,9], 1, [1,2,4], 1:14, 2:5, 2:6, 1, 1, 1, 1, 1, 4, 1});
jumping = struct( 'name', 'jumping', 'subject', {13,16,49} , 'trial', {[11,13,19,32,39:42], [1:7,9:10], 2:3});
varied = struct( 'name', 'varied', 'subject', {49} , 'trial', {4:5});
basketball = struct( 'name', 'basketball' , 'subject', {6} , 'trial', {2:15});
dance = struct( 'name', 'dance', 'subject', {5,49,55,93} , 'trial', {2:20, [9:17,22], 1:2, 6});
gymnastics = struct( 'name', 'gymnastics', 'subject', {49} , 'trial', {6:8});
acrobatics = struct( 'name', 'acrobatics', 'subject', {49} , 'trial', {21});
martial_arts = struct( 'name', 'martial_arts', 'subject', {2,12} , 'trial', {[5,7:9], 4});
soccer = struct( 'name', 'soccer', 'subject', {10,11} , 'trial', {1:6, 1});
boxing = struct( 'name', 'boxing', 'subject', {13,14,15,17} , 'trial', {17:18, 1:3, 13, 10});
general_exercise_and_stretching = struct( 'name', 'general_exercise_and_stretching', 'subject', {13,14,42,49} , 'trial', {29:31, [6,14,20], 1, 18:20});
common_behaviors_and_expressions = struct( 'name', 'common_behaviors_and_expressions', 'subject', {2,13,14,15,26,40,41} , 'trial', {[6,10], [1:10,12,14:16,20:25,33:38], [4,5,7:13,15:19,21:23,27:37], [2,6,7,10], 9:11, 6:11, 7:9});
pantomime = struct( 'name', 'pantomime', 'subject', {24,25,26,27,28,29,30,31,32,54,55,56} , 'trial', {1, 1, 3:8, [3:5,7:11], [2:4,6:19], [3:6,8:25], [2:6,8:23], [2:4,6:21], [4:6,8:22], [1:19,23:27], [3,5:28], 2:8});
communication_gestures_and_signals = struct( 'name', 'communication_gestures_and_signals', 'subject', {13,14,15,26,27,28,29,30,31,32,40,41} , 'trial', {26:28, 24:26, 8, 2, [2,6], [1,5], [2,7], [1,7], [1,5], [3,7], 12, 10:11});
cross_category_variety = struct( 'name', 'cross_category_variety', 'subject', {15} , 'trial', {[4,5,11,12]});

% Category structs
human_interaction = struct( 'subcategoryNames', {'two_subjects'}, 'subcategories', {two_subjects});
interaction_with_environment = struct( 'subcategoryNames', {'playground','uneven_terrain','path_with_obstacles'}, 'subcategories', {playground, uneven_terrain, path_with_obstacles});
locomotion = struct( 'subcategoryNames', {'running', 'walking', 'jumping', 'varied'}, 'subcategories', {running, walking, jumping, varied});
physical_activities_and_sports = struct('subcategoryNames', {'basketball', 'dance', 'gymnastics', 'acrobatics', 'martial_arts', 'soccer', 'boxing', 'general_exercise_and_stretching'}, ...
    'subcategories', {basketball, dance, gymnastics, acrobatics, martial_arts, soccer, boxing, general_exercise_and_stretching});
situations_and_scenarios = struct('subcategoryNames', {'common_behaviors_and_expressions', 'pantomime', 'communication_gestures_and_signals', 'cross_category_variety'},...
    'subcategories', {common_behaviors_and_expressions, pantomime, communication_gestures_and_signals, cross_category_variety});

% Use this struct as a argument in collectSkeletonData() to load all motions.
categories = {human_interaction, interaction_with_environment, locomotion, physical_activities_and_sports, situations_and_scenarios}; 