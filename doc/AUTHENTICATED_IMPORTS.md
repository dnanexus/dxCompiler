# Authenticated HTTP Imports for WDL

This document describes how to import WDL files from HTTP sources that require authentication.

## Overview

dxCompiler can import WDL files from URLs that require authentication, such as:
- Private GitHub repositories
- Private GitLab repositories
- Internal corporate servers

Authentication is provided via per-domain Bearer tokens set in a single environment variable.

## Environment Variable

| Variable | Required | Description |
|----------|----------|-------------|
| `WDL_IMPORT_TOKENS` | No | Semicolon-separated `domain:token` pairs for authenticated HTTP imports |

**Format:** `domain:token[;domain:token]*`

Tokens are only sent to domains explicitly listed in this variable. Requests to unlisted domains proceed without authentication.

## Configuration Examples

### Basic GitHub Private Repository Access

```bash
# Generate a token at https://github.com/settings/tokens
# Required scope: repo (for private repositories)
export WDL_IMPORT_TOKENS="raw.githubusercontent.com:ghp_xxxxxxxxxxxxxxxxxxxx"

java -jar dxCompiler.jar compile workflow.wdl -project project-xxxx -folder /my/workflows/
```

### Multiple Private Sources

```bash
# Different tokens for different services
export WDL_IMPORT_TOKENS="raw.githubusercontent.com:ghp_xxxxxxxxxxxxxxxxxxxx;gitlab.com:glpat-yyyyyyyyyy;internal.company.com:my-internal-token"

java -jar dxCompiler.jar compile workflow.wdl -project project-xxxx -folder /my/workflows/
```

## WDL Import Syntax

### GitHub Raw Content URL

```wdl
import "https://raw.githubusercontent.com/owner/repo/branch/path/to/file.wdl"
```

### GitHub Blob URL (Not Recommended)

GitHub blob URLs (`github.com/owner/repo/blob/...`) do not return raw content.
Use raw.githubusercontent.com URLs instead.

## Error Messages

### 401 Unauthorized

```
HTTP 401 Unauthorized when accessing https://raw.githubusercontent.com/...
If this is a private repository, ensure WDL_IMPORT_TOKENS is set.
Format: domain:token[;domain:token]*
```

**Solution:** Set the `WDL_IMPORT_TOKENS` environment variable with the appropriate domain and token.

### 403 Forbidden

```
HTTP 403 Forbidden when accessing https://raw.githubusercontent.com/...
The token may be invalid or lack the required permissions.
```

**Solution:** The token may be invalid or lack the required permissions. For GitHub, ensure the token has the `repo` scope.

## Security Considerations

1. **Token Scope**: Only grant the minimum required permissions to your token
2. **Per-Domain Tokens**: Each domain gets its own token; tokens are never sent to domains they aren't configured for
3. **No Logging**: Token values are never logged; only domain names are traced
4. **HTTPS Recommended**: Always use HTTPS URLs for private imports
5. **Token Format**: Tokens containing colons are supported (only the first colon in each entry is used as a delimiter). Tokens containing semicolons are not supported.

## Debugging

Enable trace logging to see when authenticated imports are used:

```bash
java -jar dxCompiler.jar compile workflow.wdl -verbose -verboseKey FileSourceResolver
```

Look for log entries like:
```
[TRACE] WDL_IMPORT_TOKENS found; authenticated HTTP imports enabled for domains: raw.githubusercontent.com, gitlab.com
[TRACE] Using authenticated HTTP for import from: raw.githubusercontent.com
```

## Troubleshooting

### Token not being sent

1. Verify `WDL_IMPORT_TOKENS` is set: `echo $WDL_IMPORT_TOKENS`
2. Check that the domain in your import URL matches a domain in the variable
3. Enable verbose logging to see authentication attempts

### Token rejected

1. Verify the token is still valid (not expired)
2. Check token permissions/scopes
3. For GitHub, ensure the token has access to the specific repository

### Public imports stopped working

The authenticated HTTP protocol is backward compatible. If `WDL_IMPORT_TOKENS` is not set, it works like the standard HTTP protocol. Verify the variable is not set if you're testing public access.

## How It Works

dxCompiler uses a custom `AuthenticatedHttpFileAccessProtocol` that:

1. Parses `WDL_IMPORT_TOKENS` into a map of domain -> token
2. For each HTTP import, looks up the domain in the map
3. If a token is found for the domain, adds an `Authorization: Bearer <token>` header to the request
4. If no token is configured for the domain, the request proceeds without authentication

This ensures that:
- Tokens are never sent to unconfigured domains
- Different services can use different tokens
- Existing workflows continue to work without modification
- Authentication failures produce clear, actionable error messages
