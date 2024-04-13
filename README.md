# Yuksel Rating
This "library" helps perform skill-based matchmaking between two parties called "players". This is based off of Cem Yuksel's (cem@cemyuksel.com) "Skill-Based Matchmaking for Competitive Two-Player Games", ACM 2024.
I am in no way associated with Cem Yuksel nor their works.

## Building
This is a single-header library, which means you `#define YUKSEL_RATING_IMPLEMENTATION` in your desired file / translation unit before including `yr.h`.
You can also make use of `yr.c` if you want your own standalone translation unit.

## Dependencies
This code depends on `math.h` and `float.h`.
No allocations required.

## License
Check the end of `yr.h`
